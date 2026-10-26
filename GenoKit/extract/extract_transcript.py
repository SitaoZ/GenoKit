# -*- coding: utf-8 -*-
'''
filename: extract_transcript.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''


import sys
import pickle
import logging
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
import multiprocessing as mp
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from GenoKit.database.database import gb_iter, genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database, transcript_count,
)
from GenoKit.utils.util import seq_upper_lower, parse_output, stream_parse_output
from GenoKit.utils.gtf_gff_tool import gtf_attr_to_gff, gff_attr_to_gtf


_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Chrom', 'Start Genome', 'End Genome',
    'Start Transcript', 'End Transcript', 'Strand','Sequence', 'Length'
]

_DNA_COMPLEMENT = str.maketrans('ATGCatgc', 'TACGtacg')
_RNA_COMPLEMENT = str.maketrans('AUGCaugc', 'UACGuacg')
_TRANSCRIPT_WRITE_BUFFER_SIZE = 1024 * 1024


'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''
def stop_codon_from_gffutils(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s

def reverse_complement(seq, is_dna=True):
    """
    sequence reverse complement
    :param seq:
    :return:
    """
    # 定义互补规则
    seq = seq.upper()
    if is_dna:
        complement = _DNA_COMPLEMENT
    else:
        complement = _RNA_COMPLEMENT

    # 反转序列并翻译互补碱基
    return seq[::-1].translate(complement)


def mature_transcript(g, t , db, genome, strand) -> str:
    '''

    :param g: gene id
    :param t: transcript id
    :param db: binary database from create
    :param genome: genoem Fasta object
    :return: mature transcript sequence
    '''
    # db' exon is coordinate by position
    seq = [genome[e.chr][e.start - 1: e.end].seq for e in db[g][t]['exon']]
    seq = ''.join(seq)
    # 5'-->3'
    if strand == '-':
        seq = reverse_complement(seq, is_dna=True)
    else:
        seq = seq.upper()


    return seq


def _feature_sequence(feature, genome, sequence_cache=None):
    """Read one feature sequence, optionally reusing it within a gene."""
    key = (feature.chr, feature.start, feature.end)
    if sequence_cache is not None:
        cached = sequence_cache.get(key)
        if cached is not None:
            return cached
    sequence = genome[feature.chr][feature.start - 1:feature.end].seq
    if sequence_cache is not None:
        sequence_cache[key] = sequence
    return sequence


def anchor_CDS(g, t, mt, db, genome, style, sequence_cache=None) -> tuple:
    '''
    parameters:
     g : gene id
     db : database from create
     transcript:  transcript id
     genome: genome Fasta object
     style: database file type
    return:
     start_codon position 
     stop_codon positiopn 
     cds length 
    '''
    mRNA_strand = db[g][t]['mRNA'].strand
    utr5, cds, utr3 = '', '', ''
    utr5_len, cds_len, utr3_len = 0, 0, 0
    # mature_transcript_seq = mature_transcript(g, t , db, genome, mRNA_strand)
    mature_transcript_seq = mt
    if style == 'gff':
        # gff : CDS contains stop codon
        if db[g][t].get('CDS'):
            cds = [_feature_sequence(c, genome, sequence_cache) for c in db[g][t]['CDS']]
            cds = ''.join(cds)
            # 5'-->3'
            cds_len = len(cds)
        if db[g][t].get('five_prime_UTR'):
            utr5 = [_feature_sequence(u, genome, sequence_cache) for u in db[g][t]['five_prime_UTR']]
            utr5_len = sum(map(len, utr5))
        if db[g][t].get('three_prime_UTR'):
            utr3 = [_feature_sequence(u, genome, sequence_cache) for u in db[g][t]['three_prime_UTR']]
            utr3_len = sum(map(len, utr3))
        # assert len(mature_transcript_seq) == utr5_len + cds_len + utr3_len
    else:
        # gtf : CDS don't contains stop codon
        if db[g][t].get('CDS'):
            # cds is coordinate by position
            cds = [_feature_sequence(c, genome, sequence_cache) for c in db[g][t]['CDS']]
            cds = ''.join(cds)
            # add stop codon
            # stop_codon_seq = ''
            # if db[g][t].get('stop_codon'):
            #     # stop codon 可能含有多个stop_codon 分布在不同的外显子上
            #     for c in db[g][t]['stop_codon']:
            #         chrom, start, end, strand = c.chr, c.start, c.end, c.strand
            #         stop_codon_seq += str(genome[chrom][start-1:end])
            stop_codon_seq = [_feature_sequence(c, genome, sequence_cache) for c in db[g][t]['stop_codon']]
            stop_codon_seq = ''.join(stop_codon_seq)
            # strand
            if mRNA_strand == '+':
                cds = cds + stop_codon_seq
            else:
                cds = stop_codon_seq + cds

            # 5'-->3'
            if mRNA_strand == '-':
                cds = reverse_complement(cds, is_dna=True)
            cds = str(cds.upper())
            cds_len = len(cds) # CDS don't contains stop codon
            utr5 = mature_transcript_seq.split(cds)[0]
            utr3 = mature_transcript_seq.split(cds)[1]
            utr5_len = len(utr5)
            utr3_len = len(utr3)

    # mt, cds2, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, style)
    # assert cds == cds2

    return utr5_len, utr3_len, cds_len





def _transcript_output_row(g, symbol, t, db, genome, style, upper, sequence_cache,
                           output_format):
    """Build one CSV row or FASTA record without temporary SeqRecord objects."""
    transcript = db[g][t]
    mrna = transcript['mRNA']
    sequence = ''.join(_feature_sequence(exon, genome, sequence_cache)
                       for exon in transcript['exon'])
    if mrna.strand == '-':
        sequence = reverse_complement(sequence, is_dna=True)

    if transcript.get('CDS'):
        utr5_len, _, cds_len = anchor_CDS(
            g, t, sequence, db, genome, style, sequence_cache)
        cds_start_transcript = utr5_len + 1
        cds_end_transcript = utr5_len + cds_len
    else:
        utr5_len = cds_len = cds_start_transcript = cds_end_transcript = 0

    if upper:
        sequence = seq_upper_lower(sequence, utr5_len, utr5_len + cds_len)
    if output_format == 'fasta':
        # Match SeqIO's titles and 60-column wrapping without per-record writers.
        identifier = t.replace('\n', ' ').replace('\r', ' ')
        description = (f'strand:{mrna.strand} start:{mrna.start} end:{mrna.end} '
                       f'length={len(sequence)} '
                       f'CDS={cds_start_transcript}-{cds_end_transcript}')
        title = (description if description.split(None, 1)[0] == identifier
                 else f'{identifier} {description}')
        return f'>{title}\n' + ''.join(
            sequence[i:i + 60] + '\n' for i in range(0, len(sequence), 60))
    return (f'{g},{symbol},{t},{mrna.chr},{mrna.start},{mrna.end},'
            f'{cds_start_transcript},{cds_end_transcript},{mrna.strand},'
            f'{sequence},{len(sequence)}\n')


def _write_all_transcripts(db, genome, output_path, style, upper, output_format):
    """Stream CSV/FASTA with gene-local sequence caching and bounded writes."""
    total = transcript_count(db)
    buffered_rows = []
    buffered_size = 0
    with open(output_path, 'w') as out:
        if output_format == 'csv':
            out.write(','.join(_CSV_HEADER) + '\n')
        with tqdm(total=total, ncols=80, desc='Transcript processing:',
                  colour='GREEN') as progress:
            current_gene = None
            sequence_cache = {}
            for gene_id, symbol, transcript_id, transcript in iter_transcript_models_with_symbols(
                db, include_feature_attributes=False
            ):
                if gene_id != current_gene:
                    sequence_cache = {}
                    current_gene = gene_id
                transcript_db = {gene_id: {'symbol': symbol, transcript_id: transcript}}
                row = _transcript_output_row(
                    gene_id, symbol, transcript_id, transcript_db, genome, style,
                    upper, sequence_cache, output_format)
                buffered_rows.append(row)
                buffered_size += len(row)
                progress.update()
                if buffered_size >= _TRANSCRIPT_WRITE_BUFFER_SIZE:
                    out.write(''.join(buffered_rows))
                    buffered_rows.clear()
                    buffered_size = 0
        if buffered_rows:
            out.write(''.join(buffered_rows))




def sub_transcript(params) -> list:
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    g, symbol, t, db, genome, style, upper, output_format = params
    # transcript_seq_list = deque()
    transcript_seq_list = []
    mRNA = db[g][t]['mRNA']
    # exon is coordinate by position
    seq = [genome[e.chr][e.start - 1: e.end].seq for e in db[g][t]['exon']]
    seq = ''.join(seq)
    seq = Seq(seq)
    if mRNA.strand == '-':
        seq = seq.reverse_complement()
    seq = seq.upper()
    # position in transcript 1 based
    if db[g][t].get('CDS'):
        utr5_len, utr3_len, cds_len = anchor_CDS(g, t, seq, db, genome, style)
        cds_start_transcript = utr5_len + 1
        cds_end_transcript = utr5_len + cds_len
    else:
        utr5_len, cds_len, cds_start_transcript, cds_end_transcript = 0, 0, 0, 0


    if upper:
        seq = seq_upper_lower(seq, utr5_len , utr5_len + cds_len)

    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand, mRNA.start,
                                                              mRNA.end, len(seq),
                                                              cds_start_transcript,
                                                              cds_end_transcript)
        transcript_record = SeqRecord(seq, id=t, description=desc)
        transcript_seq_list.append(transcript_record)
    else:
        seq = str(seq) # for polars write
        it = [g, symbol, t, mRNA.chr, mRNA.start, mRNA.end,
              cds_start_transcript, cds_end_transcript, mRNA.strand, seq, len(seq)]
        transcript_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return transcript_seq_list


def sub_transcript_faster(params) -> list:
    """
        parameters:
            params: argparse
        return:
            transcript_seq_list, a list
        """
    g, symbol, t, db, genome, style, upper, output_format = params
    transcript_seq_list = []
    mRNA = db[g][t]['mRNA']
    # print("mRNA:", mRNA) # 使用msgpack会导致namedtupled失效,转变成普通的tuple
    mRNA_strand = mRNA.strand
    # mt: exon is coordinate by position
    seq = [genome[e.chr][e.start - 1: e.end].seq for e in db[g][t]['exon']]
    seq = ''.join(seq)
    if mRNA_strand == '-':
        seq = reverse_complement(seq, is_dna=True)
    # position in transcript 1 based
    if db[g][t].get('CDS'):
        utr5_len, utr3_len, cds_len = anchor_CDS(g, t, seq, db, genome, style)
        cds_start_transcript = utr5_len + 1
        cds_end_transcript = utr5_len + cds_len
    else:
        utr5_len, cds_len, cds_start_transcript, cds_end_transcript = 0, 0, 0, 0
    if upper:
        seq = seq_upper_lower(seq, utr5_len, utr5_len + cds_len)

    if output_format == 'fasta':
        desc = 'strand:%s start:%d end:%d length=%d CDS=%d-%d' % (mRNA_strand,
                                                                  mRNA.start,
                                                                  mRNA.end,
                                                                  len(seq),
                                                                  cds_start_transcript,
                                                                  cds_end_transcript)

        seq = Seq(seq)
        transcript_record = SeqRecord(seq, id=t, description=desc)
        transcript_seq_list.append(transcript_record)

    elif output_format == 'gff':
        if "=" not in mRNA.attribute:
            # gtf to gff
            # attribute = gtf_attr_to_gff(mRNA.attribute)
            attribute = f"ID={t};Gene={g};"
            attribute = attribute.replace("gene_name", "Name")
        else:
            attribute = mRNA.attribute
        gff_format = f"{mRNA.chr}\t{mRNA.source}\t{mRNA.type}\t{mRNA.start}\t{mRNA.end}\t{mRNA.score}\t{mRNA.strand}\t{mRNA.phase}\t{attribute}"
        transcript_seq_list.append(gff_format)

    elif output_format == 'gtf':
        if "=" in mRNA.attribute:
            # gff to gtf
            # attribute = gff_attr_to_gtf(mRNA.attribute)
            attribute = f'gene_id {g}; transcript_id {t};'
        else:
            attribute = mRNA.attribute

        gtf_format = f"{mRNA.chr}\t{mRNA.source}\t{mRNA.type}\t{mRNA.start}\t{mRNA.end}\t{mRNA.score}\t{mRNA.strand}\t{mRNA.phase}\t{attribute}"
        transcript_seq_list.append(gtf_format)

    else:
        transcript_seq_list.append([g, symbol, t, mRNA.chr, mRNA.start, mRNA.end,
                                    cds_start_transcript, cds_end_transcript,
                                    mRNA_strand, seq, len(seq)])
    return transcript_seq_list
    
def get_transcript(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Transcript.log"
    )
    logging.info("Loading database")
    import time
    # starttime = time.time()
    # method1
    db, t2g = load_annotation_database(args.database)

    # method2
    # import larch.pickle as lp
    # db, t2g = lp.load(f)

    # method3
    # import msgpack
    # db, t2g = msgpack.load(f)

    # method4
    # file_style = args.style.lower()
    # from GenoKit.utils.util import gff_feature_dict
    # from GenoKit.utils.util import gtf_feature_dict
    # if file_style == 'gff':
    #     db, t2g = gff_feature_dict(args.database, file_style)
    # elif file_style == 'gtf':
    #     db, t2g = gtf_feature_dict(args.database, file_style)

    # endtime = time.time()
    # print("Time:", endtime - starttime)

    # Fasta object
    logging.info("Loading genome")
    genome = genome_dict(args.genome)

    if not args.transcript:
        if args.output and (args.output_format == 'csv' or
                            (args.output_format == 'fasta' and not args.print)):
            logging.info("Transcript processing at whole-genome scale")
            _write_all_transcripts(
                db, genome, args.output, args.style, args.upper, args.output_format)
            logging.info("Transcript finished at whole-genome scale")
            return

        logging.info("Transcript processing at whole-genome scale")
        param_iter = (
            (
                g, gene.get('symbol') or get_gene_symbol(db, g), t, db,
                genome, args.style, args.upper, args.output_format,
            )
            for g, gene in db.items()
            for t in gene if t not in {'gene', 'symbol'}
        )
        with open(args.output, 'w') as out:
            if args.output_format == 'csv':
                out.write(",".join(_CSV_HEADER) + '\n')
            for para in tqdm(param_iter, ncols = 80, total=transcript_count(db, t2g), \
                             desc='Transcript processing:', colour="GREEN"):
                transcript_seq_list = sub_transcript_faster(para)
                # if transcript_seq_list == None:
                #     continue
                stream_parse_output(args, transcript_seq_list, out)
        logging.info("Transcript finished at whole-genome scale")
    else:
        logging.info("Transcript processing at transcript scale")
        t = args.transcript
        if not t2g.get(t):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        g = t2g[t]
        param = (
            g, get_gene_symbol(db, g), t, db, genome, args.style,
            args.upper, args.output_format,
        )
        transcript_seq_list = sub_transcript(param)
        parse_output(args, transcript_seq_list)
        logging.info("Transcript finished at transcript scale")


def get_transcript_parallel(args):
    '''
    multiple processes
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db, t2g = load_annotation_database(args.database)
    # Fasta object，不能作为多进程的传递的参数
    genome = genome_dict(args.genome)
    if not args.transcript:
        param_list = (
            (
                g, gene.get('symbol') or get_gene_symbol(db, g), t, db,
                args.genome, args.style, args.upper, args.output_format,
            )
            for g, gene in db.items()
            for t in gene if t not in {'gene', 'symbol'}
        )
        # transcript_seq_list = deque()
        with mp.Pool(processes=args.process) as pool:
            # 使用tqdm显示进度
            transcript_seq_list = list(tqdm(pool.imap(sub_transcript, param_list),
                                total=transcript_count(db, t2g),
                                ncols=100,
                                desc=f'Transcript processing ({args.process} cpu mode):'
                                ))
        transcript_seq_list = [d for de in transcript_seq_list if de != None for d in de]
        parse_output(args, transcript_seq_list)
    else:
        t = args.transcript
        g = t2g[t]
        param = (
            g, get_gene_symbol(db, g), t, db, genome, args.style,
            args.upper, args.output_format,
        )
        transcript_seq_list = sub_transcript_faster(param)
        parse_output(args, transcript_seq_list)


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_transcript_gb(args):
    '''
    '''
    transcript = []
    for record in gb_iter(args.genbank):
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                transcript_seq = ''
                #for part in feature.location.parts:
                #    cds_seq += part.extract(record.seq)
                transcript_seq = feature.location.extract(record).seq
                # part.strand 会将FeatureLocation -1的反向互补
                # geneid or locus_tag 
                if 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                else:
                    gene_id = "Null"
                # protein id  
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                else:
                    protein_id = 'Null'

                transcript_seq_record = SeqRecord(transcript_seq,id='gene:%s protein:%s'%(gene_id, protein_id),
                                 description='strand %s length %d'%(feature.strand, len(transcript_seq)))
                transcript.append(transcript_seq_record)
                break 
        break     
    if args.print and args.format == 'dna':
        SeqIO.write(transcript, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(transcript, sys.stdout, 'fasta')
    elif args.output and args.format == 'dna':
        SeqIO.write(transcript, args.output, 'fasta')
    elif args.output and args.format == 'protein':
        SeqIO.write(transcript, args.output, 'fasta')
