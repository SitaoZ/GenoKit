#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: kmer_stat.py
date: 2025/12/1 09:42
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from tqdm import tqdm
import _pickle as cPickle
import multiprocessing as mp
from Bio.Seq import Seq
from pyfaidx import Fasta, Faidx
from collections import defaultdict, Counter
from GenoKit.database.fasta_split import fasta_count
from GenoKit.database.kmer_database import (
    JellyfishKmerCounts,
    register_jellyfish_dataset,
    register_sqlite_counts,
)


_PYTHON_MAX_FASTA_BYTES = 250 * 1024 * 1024
_BASE_TO_BITS = [-1] * 256
for _base, _bits in ((ord('A'), 0), (ord('C'), 1), (ord('G'), 2), (ord('T'), 3),
                     (ord('a'), 0), (ord('c'), 1), (ord('g'), 2), (ord('t'), 3)):
    _BASE_TO_BITS[_base] = _bits
_BITS_TO_BASE = 'ACGT'


def _kmer_output_path(args, suffix):
    basename = getattr(args, '_output_basename', None) or os.path.basename(args.fasta)
    return os.path.join(args.output_dir, f'{basename}.kmer.{args.length}.{suffix}')


def _resolve_jellyfish(args):
    if args.jellyfish:
        candidate = os.path.expanduser(args.jellyfish)
        if os.path.sep in candidate:
            return candidate if os.path.isfile(candidate) and os.access(candidate, os.X_OK) else None
        return shutil.which(candidate)
    return shutil.which('jellyfish')


def _run_jellyfish(args, executable):
    """Count k-mers with Jellyfish, the scalable path for genome-sized FASTA."""
    output_path = _kmer_output_path(args, 'jf')
    threads = args.threads or (os.cpu_count() or 1)
    command = [
        executable, 'count',
        '-m', str(args.length),
        '-s', str(args.hash_size),
        '-t', str(threads),
        '-o', output_path,
    ]
    if args.reverse_complement:
        command.append('-C')
    command.append(args.fasta)
    print(f'Counting k-mers with Jellyfish ({threads} threads)...')
    subprocess.run(command, check=True)
    print(f'Jellyfish database written to: {output_path}')
    register_jellyfish_dataset(
        args.sqlite_db,
        output_path,
        args.length,
        args.source_type,
        args.reverse_complement,
        getattr(args, '_source_description', args.fasta),
        executable,
    )
    print(
        f'K-mer dataset registered in: {args.sqlite_db} '
        f'(source={args.source_type}, k={args.length})'
    )
    if args.reverse_complement:
        print('Reverse-complement k-mers were counted in canonical mode.')
    return output_path


def _decode_kmer(value, length):
    characters = [''] * length
    for index in range(length - 1, -1, -1):
        characters[index] = _BITS_TO_BASE[value & 3]
        value >>= 2
    return ''.join(characters)


def _count_kmers_python(args):
    """Exact, streaming fallback for small FASTA files.

    It uses a rolling two-bit representation rather than allocating a string
    slice for every position.  Ambiguous-base windows are skipped, which is the
    standard behavior for genome k-mer counting.
    """
    fasta_size = os.path.getsize(args.fasta)
    if fasta_size > _PYTHON_MAX_FASTA_BYTES:
        size_mb = fasta_size / (1024 * 1024)
        raise SystemExit(
            f'Python k-mer counting is limited to FASTA files up to '
            f'{_PYTHON_MAX_FASTA_BYTES // (1024 * 1024)} MB (input: {size_mb:.0f} MB). '
            'Install Jellyfish and run with --engine jellyfish or -j PATH for genome-scale counting.'
        )
    if not 1 <= args.length <= 31:
        raise SystemExit('The Python engine supports k-mer lengths from 1 to 31.')

    mask = (1 << (2 * args.length)) - 1
    reverse_shift = 2 * (args.length - 1)
    counts = Counter()
    forward = reverse = valid_length = 0
    with open(args.fasta, 'rb') as handle, tqdm(
        total=fasta_size,
        unit='B',
        unit_scale=True,
        mininterval=0.5,
        ncols=80,
        desc='Counting k-mers',
        colour='GREEN',
    ) as progress:
        for line in handle:
            progress.update(len(line))
            if line.startswith(b'>'):
                forward = reverse = valid_length = 0
                continue
            for character in line.rstrip(b'\r\n'):
                bits = _BASE_TO_BITS[character]
                if bits < 0:
                    forward = reverse = valid_length = 0
                    continue
                forward = ((forward << 2) | bits) & mask
                reverse = (reverse >> 2) | ((3 - bits) << reverse_shift)
                valid_length += 1
                if valid_length < args.length:
                    continue
                if args.reverse_complement:
                    counts[min(forward, reverse)] += 1
                else:
                    counts[forward] += 1

    output_path = _kmer_output_path(args, 'bin')
    # Preserve the legacy pickle-dictionary format for downstream callers.
    decoded_counts = {_decode_kmer(kmer, args.length): count for kmer, count in counts.items()}
    with open(output_path, 'wb') as handle:
        cPickle.dump(decoded_counts, handle, protocol=-1)
    register_sqlite_counts(
        args.sqlite_db,
        decoded_counts,
        args.length,
        args.source_type,
        args.reverse_complement,
        getattr(args, '_source_description', args.fasta),
    )
    print(f'Counted {len(decoded_counts):,} unique k-mers.')
    print(f'K-mer dictionary written to: {output_path}')
    print(
        f'K-mer counts stored in: {args.sqlite_db} '
        f'(source={args.source_type}, k={args.length})'
    )
    return decoded_counts





# multiple processes
def process_chromosome(args):
    """
    处理单个染色体的kmer统计
    """
    chrom_name, fasta_path, kmer_length, rc = args
    kmer_bin = defaultdict(int)

    try:
        # genome file format only support fasta an BGZIP file.
        genome = Fasta(fasta_path, rebuild=False, build_index=False)
        sequence = genome[chrom_name]

        if rc:
            for i in range(0, len(sequence), 1):
                if i + kmer_length > len(sequence):
                    continue
                kmer = sequence[i : i+kmer_length].seq.upper()
                # https://github.com/mdshw5/pyfaidx
                rc_kmer = sequence[i : i+kmer_length].reverse.complement.seq.upper()
                kmer_bin[kmer] += 1
                kmer_bin[rc_kmer] += 1
        else:
            for i in range(0, len(sequence), 1):
                if i + kmer_length > len(sequence):
                    continue
                kmer = sequence[i : i+kmer_length].seq.upper()
                # 过滤掉包含非ATCGN的kmer
                #if all(base in 'ATCGN' for base in kmer):
                kmer_bin[kmer] += 1

    except Exception as e:
        print(f"处理染色体 {chrom_name} 时出错: {e}")

    return dict(kmer_bin)


def merge_results(results):
    """
    合并所有进程的结果
    """
    merged = defaultdict(int)
    for result in results:
        for kmer, count in result.items():
            merged[kmer] += count
    return dict(merged)


def get_kmer_parallel_no_tdpm(args):
    """
    使用多进程进行kmer分析
    """
    # 构建索引
    if not os.path.exists(args.fasta + ".fai"):
        # print("构建FASTA索引...")
        fa = Faidx(args.fasta)
        fa.build_index()

    # 获取所有染色体名称
    # genome file format only support fasta an BGZIP file.
    genome = Fasta(args.fasta, rebuild=False, build_index=False)
    chrom_names = list(genome.keys())
    # print(f"检测到 {len(chrom_names)} 条染色体/contigs")

    # 准备进程参数
    process_args = [(chrom, args.fasta, args.length) for chrom in chrom_names]

    # 使用多进程处理
    # print(f"使用 {args.processes} 个进程进行处理...")
    with mp.Pool(processes=args.processes) as pool:
        results = pool.map(process_chromosome, process_args)

    # 合并结果
    # print("合并结果...")
    kmer_bin = merge_results(results)

    # 保存结果
    kmer_bin_path = os.path.join(args.output_dir, f'{os.path.basename(args.fasta)}.kmer.{args.length}.bin')
    with open(kmer_bin_path, 'wb') as f:
        cPickle.dump(kmer_bin, f, protocol=-1)

    # print(f"kmer统计完成! 共找到 {len(kmer_bin)} 个独特的kmer")
    return kmer_bin


def get_kmer_parallel(args):
    """
    带进度显示的多进程版本 tqdm
    """
    if not os.path.exists(args.fasta + ".fai"):
        # print("构建FASTA索引...")
        fa = Faidx(args.fasta)
        fa.build_index()

    genome = Fasta(args.fasta, rebuild=False, build_index=False)
    chrom_names = list(genome.keys())
    # 对基因组进行kmer分析 或者 对基因的成熟转录本进行kmer分析
    process_args = [(chrom, args.fasta, args.length, args.reverse_complement) for chrom in chrom_names]
    # print(f"使用 {args.process} 个进程处理 {len(chrom_names)} 条染色体...")
    # args.process
    with mp.Pool(processes=2) as pool:
        # 使用tqdm显示进度
        results = list(tqdm(pool.imap(process_chromosome, process_args),
                            total=len(chrom_names),
                            ncols=80,
                            desc=f'Kmer processing ({2} cpu mode):'
                            ))

    kmer_bin = merge_results(results)

    kmer_bin_path = os.path.join(args.output_dir, f'{os.path.basename(args.fasta)}.kmer.{args.length}.bin')
    with open(kmer_bin_path, 'wb') as f:
        cPickle.dump(kmer_bin, f, protocol=-1)

    # print(f"kmer统计完成! 共找到 {len(kmer_bin)} 个独特的kmer")
    return kmer_bin



# single process
def get_kmer_single(args):
    """
    单进程版本
    """
    if not os.path.exists(args.fasta + ".fai"):
        fa = Faidx(args.fasta)
        fa.build_index()

    # kmer_bin = defaultdict(int)
    kmer_bin = Counter()
    genome = Fasta(args.fasta, rebuild=False, build_index=False, strict_bounds=False)
    chrom_names = list(genome.keys())
    for chrom in tqdm(genome.keys(), ncols = 80, total=len(chrom_names),
                      desc='Kmer processing:', colour="GREEN"):

        sequence = genome[chrom]
        if args.reverse_complement:
            for i in range(0, len(sequence), 1):
                # if i + args.length > len(sequence):
                #     continue
                kmer = sequence[i: i + args.length].seq.upper()
                rc_kmer = kmer.reverse.complement
                kmer_bin[kmer] += 1
                kmer_bin[rc_kmer] += 1
        else:
            for i in range(0, len(sequence), 1):
                kmer = sequence[i : i+args.length].seq.upper()
                # remove none atcgn base
                #if all(base in 'ATCGN' for base in kmer):
                kmer_bin[kmer] += 1

    kmer_bin_path = os.path.join(args.output_dir, f'{os.path.basename(args.fasta)}.kmer.{args.length}.bin')
    with open(kmer_bin_path, 'wb') as f:
        cPickle.dump(dict(kmer_bin), f, protocol=-1)

    # print(f"kmer统计完成! 共找到 {len(kmer_bin)} 个独特的kmer")
    return kmer_bin

# jellyfish
def get_kmer_jellyfish(args):
    kmer_length = args.length
    fasta_path = args.fasta
    output_dir = args.output_dir
    basename = os.path.basename(fasta_path)

    if args.jellyfish:
        cmd = f"jellyfish count -m {kmer_length} -s 3G -t 8 -o {output_dir}/{basename}.kmer.{kmer_length}.bin {fasta_path}"
        os.system(cmd)


def get_kmer_khmer(args):
    """

    :param args:
    :return:
    """
    kmer_length = args.length
    fasta_path = args.fasta
    output_dir = args.output_dir
    basename = os.path.basename(fasta_path)
    try:
        import khmer
        raise ValueError("khmer not exist, please pip install khmer.")
    except Exception as e:
        print(e)

def get_kmer_oxli(args):
    """
    kmer analysis using oxli (https://pypi.org/project/oxli/)
    :param args:
    :return:
    """
    import screed
    from oxli import KmerCountTable

    kmer_length = args.length
    fasta_path = args.fasta
    output_dir = args.output_dir
    basename = os.path.basename(fasta_path)
    database_path = f"{output_dir}/{basename}.kmer.{kmer_length}.bin"
    line_count = fasta_count(fasta_path)
    counts = KmerCountTable(ksize=kmer_length)
    with screed.open(fasta_path) as seqfile:
        for record in tqdm(seqfile, ncols = 80, total=line_count, \
                           desc="Kmer processing:", colour="GREEN"):
            # print(record .name, record .sequence)
            counts.consume(record.sequence)
    # print(counts.get("AGAGTCATTGGCAGCTGT"))
    counts.save(database_path)

def get_kmer_danio(args):
    import dnaio
    from collections import defaultdict
    k = args.length
    kmer_bin = defaultdict(int)
    with dnaio.open(args.fasta) as reader:
        for record in reader:
            seq = record.sequence
            for i in range(len(seq) - k + 1):
                kmer_bin[seq[i:i + k]] += 1
        return kmer_bin


    kmer_bin_path = os.path.join(args.output_dir, f'{os.path.basename(args.fasta)}.kmer.{args.length}.bin')
    with open(kmer_bin_path, 'wb') as f:
        cPickle.dump(dict(kmer_bin), f, protocol=-1)

    # print(f"kmer统计完成! 共找到 {len(kmer_bin)} 个独特的kmer")
    return kmer_bin

def _parse_lengths(args):
    values = []
    text = getattr(args, 'lengths', None)
    if text:
        for item in text.split(','):
            item = item.strip()
            if not item:
                continue
            if '-' in item:
                lower, upper = item.split('-', 1)
                lower, upper = int(lower), int(upper)
                if lower > upper:
                    raise SystemExit(f'Invalid k-mer range: {item}')
                values.extend(range(lower, upper + 1))
            else:
                values.append(int(item))
    else:
        values.append(int(args.length))
    values = sorted(set(values))
    if not values or any(value < 1 for value in values):
        raise SystemExit('All k-mer lengths must be positive integers.')
    return values


def _write_mature_transcript_fasta(database_path, genome_path, output_path):
    """Write every annotated mature transcript once in transcript orientation."""
    from GenoKit.database.annotation_db import load_annotation_database
    from GenoKit.database.database import genome_dict

    database, _t2g = load_annotation_database(database_path)
    genome = genome_dict(genome_path)
    transcript_count = 0
    with open(output_path, 'w', encoding='utf-8') as output:
        for gene_id in tqdm(
            database,
            total=len(database),
            ncols=80,
            desc='Building transcriptome',
            colour='GREEN',
        ):
            gene_model = database[gene_id]
            for transcript_id in gene_model:
                if transcript_id in {'gene', 'symbol'}:
                    continue
                transcript = gene_model[transcript_id]
                if 'mRNA' not in transcript or not transcript.get('exon'):
                    continue
                exons = sorted(transcript['exon'], key=lambda feature: (feature.start, feature.end))
                sequence = ''.join(
                    str(genome[exon.chr][exon.start - 1:exon.end]) for exon in exons
                )
                if transcript['mRNA'].strand == '-':
                    sequence = str(Seq(sequence).reverse_complement())
                sequence = sequence.upper()
                if not sequence:
                    continue
                output.write(f'>{transcript_id} gene={gene_id}\n')
                for start in range(0, len(sequence), 80):
                    output.write(sequence[start:start + 80] + '\n')
                transcript_count += 1
    return transcript_count


@contextmanager
def _kmer_input(args):
    """Resolve a genome or mature-transcript FASTA for the counting engines."""
    if args.source_type == 'genome':
        if not args.fasta:
            raise SystemExit('--fasta is required when --source-type genome is used.')
        args._source_description = os.path.abspath(args.fasta)
        yield args.fasta
        return

    if args.fasta:
        args._source_description = os.path.abspath(args.fasta)
        yield args.fasta
        return

    if not args.database or not args.genome:
        raise SystemExit(
            'Transcriptome k-mers require either a mature-transcript --fasta, or both '
            '--database and --genome so GenoKit can construct mature transcripts.'
        )
    temporary = tempfile.NamedTemporaryFile(
        prefix='genokit_mature_transcripts_', suffix='.fa',
        dir=args.output_dir, delete=False,
    )
    temporary_path = temporary.name
    temporary.close()
    try:
        count = _write_mature_transcript_fasta(args.database, args.genome, temporary_path)
        if count == 0:
            raise SystemExit('No mature transcript sequences could be built from the annotation database.')
        print(f'Built {count:,} mature transcript sequences for k-mer counting.')
        args._source_description = (
            f'database={os.path.abspath(args.database)};genome={os.path.abspath(args.genome)}'
        )
        args._output_basename = 'mature_transcripts'
        yield temporary_path
    finally:
        if hasattr(args, '_output_basename'):
            del args._output_basename
        try:
            os.unlink(temporary_path)
        except FileNotFoundError:
            pass


def get_kmer(args):
    """Create source-aware k-mer datasets and register them in one SQLite catalogue."""
    os.makedirs(args.output_dir, exist_ok=True)
    args.sqlite_db = os.path.abspath(
        os.path.expanduser(
            args.sqlite_db or os.path.join(args.output_dir, 'genokit.kmer.sqlite')
        )
    )
    lengths = _parse_lengths(args)
    executable = _resolve_jellyfish(args)
    if getattr(args, 'register_jf', None):
        if not executable:
            raise SystemExit(
                'Jellyfish is required to inspect existing databases; pass --jellyfish PATH.'
            )
        if not args.reverse_complement:
            raise SystemExit(
                'Use -r when registering existing Jellyfish files to confirm that they '
                'were built with canonical reverse-complement counting.'
            )
        registered = []
        for payload in args.register_jf:
            adapter = JellyfishKmerCounts(payload, jellyfish=executable)
            adapter.close()
            register_jellyfish_dataset(
                args.sqlite_db,
                payload,
                adapter.k,
                args.source_type,
                True,
                os.path.abspath(payload),
                executable,
            )
            registered.append(payload)
            print(
                f'Registered {payload} in {args.sqlite_db} '
                f'(source={args.source_type}, k={adapter.k})'
            )
        return registered
    original_fasta = args.fasta
    original_length = args.length
    results = []
    with _kmer_input(args) as fasta_path:
        args.fasta = fasta_path
        try:
            for length in lengths:
                args.length = length
                if args.engine == 'jellyfish':
                    if not executable:
                        raise SystemExit(
                            'Jellyfish was not found. Install it or pass its path with --jellyfish PATH.'
                        )
                    results.append(_run_jellyfish(args, executable))
                elif args.engine == 'auto' and executable:
                    results.append(_run_jellyfish(args, executable))
                else:
                    if args.engine == 'auto':
                        print(
                            'Jellyfish was not found; using the exact Python fallback for a small FASTA.',
                            file=sys.stderr,
                        )
                    results.append(_count_kmers_python(args))
        finally:
            args.fasta = original_fasta
            args.length = original_length
    return results[0] if len(results) == 1 else results
