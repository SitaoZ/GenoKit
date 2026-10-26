# -*- coding: utf-8 -*-
'''
filename: extract_cds.py
date: 2025/08/13 09:27
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.utils.util import stop_codon_seq
from GenoKit.utils.util import add_stop_codon
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database,
)
from GenoKit.utils.util import parse_output, stream_parse_output
from GenoKit.utils.gtf_gff_tool import gff_attr_to_gtf, gtf_attr_to_gff


'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''

# 'Genome Position', 'CDS Split'
_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Chrom',
    'Start', 'End', 'Strand', 'CDS', 'CDS Length'
]

def sub_cds(params):
    if len(params) == 6:
        g, t, db, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, db, genome, style, output_format = params
    cds_seq_list = []
    seq = ''
    cds_start_transcript = 0
    cds_end_transcript = 0
    # cds_gff_lines = []
    mRNA = db[g][t]['mRNA']
    cds_features = db[g][t].get('CDS', ())
    if cds_features:
        if mRNA.strand == '+':
            cds_start_transcript = cds_features[0].start - mRNA.start
        else:
            cds_start_transcript = mRNA.end - cds_features[-1].end
        for c in cds_features:
            # cds_gff_lines.append(c)
            chrom, start, end, strand = c.chr, c.start, c.end, c.strand
            s = str(genome[chrom][start-1:end])
            seq += s
            cds_end_transcript += len(s)
        cds_end_transcript = cds_end_transcript + cds_start_transcript
    else:
        return
    # deal with start and stop codon in gtf file, but start codon is included in CDS in GTF
    if style.lower() == 'gtf':
        # add stop codon
        stop_codon_sequence = stop_codon_seq(db, g, t, genome[chrom])
        seq = add_stop_codon(seq, mRNA.strand, stop_codon_sequence)

    seq = Seq(seq)
    if mRNA.strand == '-':
        seq = seq.reverse_complement()

    if output_format == 'gff':
        cds_gff_lines = (c for c in cds_features)
        for line in cds_gff_lines:
            array = list(line)
            if "=" not in array[8]:  # gtf to gff
                array[8] = gtf_attr_to_gff(array[8])
            # gff to gff no transfer
            gff_line = "\t".join(map(str, array))
            cds_seq_list.append(gff_line)
    elif output_format == 'gtf':
        # 需要修改
        cds_gtf_lines = (c for c in cds_features)
        for line in cds_gtf_lines:
            array = list(line)
            if "=" in array[8]:  # gff to gtf
                array[8] = gff_attr_to_gtf(array[8])
            # gtf to gtf no transfer
            gtf_line = "\t".join(map(str, array))
            cds_seq_list.append(gtf_line)

    elif output_format == 'fasta':
        # fasta (defalut)
        desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand,
                                                              mRNA.start,
                                                              mRNA.end,
                                                              len(seq),
                                                              cds_start_transcript,
                                                              cds_end_transcript)
        cdsRecord = SeqRecord(seq, id=t.replace('transcript:',''), description=desc)
        cds_seq_list.append(cdsRecord)
    else:
        # csv
        it = [g, t, mRNA.chr, mRNA.start, mRNA.end, mRNA.strand, seq, len(seq)]
        header = _CSV_HEADER
        if symbol is not None:
            it.insert(1, symbol)
        else:
            header = [column for column in _CSV_HEADER if column != 'Gene Symbol']
        cds_seq_list.append(dict(zip(header, it)))
    return cds_seq_list


def sub_cds_faster(params):
    """
    extract cds
    :param params:
    g (str): gene
    t (str): transcript
    db (dict):database from create
    genome: Fasta object
    style (str): gtf/gff
    output_format (str): csv, gff, fasta
    :return:
    cds_seq_list List[Dict]
    """
    if len(params) == 6:
        g, t, database_or_transcript, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, database_or_transcript, genome, style, output_format = params
    if 'mRNA' in database_or_transcript:
        transcript = database_or_transcript
    else:
        transcript = database_or_transcript[g][t]
    cds_seq_list = []
    cds_features = transcript['CDS']
    if output_format == 'gff':
        cds_gff_lines = (c for c in cds_features)
        for cds_index, cds in enumerate(cds_gff_lines, 1):
            if cds.strand == '-':
                cds_index = len(cds_features) - cds_index + 1
            array = list(cds)
            # Jbrowser2 needs Parent
            array[8] = f'ID={t}:cds:{cds_index};Parent={t};GeneID={g};'
            gff_line = "\t".join(map(str, array))
            cds_seq_list.append(gff_line)
    elif output_format == 'gtf':
        cds_gtf_lines = (c for c in cds_features)
        for cds_index,cds in enumerate(cds_gtf_lines, 1):
            if cds.strand == '-':
                cds_index = len(cds_features) - cds_index + 1
            array = list(cds)
            array[8] = f'gene_id "{g}"; transcript_id "{t}"; cds_id "{t}:cds:{cds_index}";'
            gtf_line = "\t".join(map(str, array))
            cds_seq_list.append(gtf_line)

    else:
        mRNA = transcript['mRNA']
        if mRNA.strand == '+':
            cds_start_position = cds_features[0].start - mRNA.start
        else:
            cds_start_position = mRNA.end - cds_features[-1].end
        seq = [genome[c.chr][c.start-1 : c.end].seq for c in cds_features]
        seq = ''.join(seq)
        cds_end_position = len(seq) + cds_start_position
        # cds genome position list
        # cds_genome_position_str = "|".join([f"{c.start}:{c.end}" for c in db[g][t]['CDS']])
        # cds_split = "|".join([genome[c.chr][c.start-1 : c.end].seq for c in db[g][t]['CDS']])

        chrom = mRNA.chr
        if not transcript['exon']:
            print(f"Warning: Gene {g} transcript {t} has no exons.")
        # deal with start and stop codon in gtf file
        # CDS contains start cond in GTF
        if style.lower() == 'gtf':
            # add stop codon
            stop_codon_sequence = stop_codon_seq(
                {g: {t: transcript}}, g, t, genome[chrom]
            )
            seq = add_stop_codon(seq, mRNA.strand, stop_codon_sequence)

        seq = Seq(seq)
        if mRNA.strand == '-':
            seq = seq.reverse_complement()
        seq = seq.upper()
        # check
        if len(seq) % 3 != 0:
            warning = f"The cds length of {t} ({len(seq)}) is not an integer multiple of 3. You need to check the annotation file."
            logging.warning(warning)

        if output_format == 'fasta':
            # fasta (defalut)
            desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand,
                                                                  mRNA.start,
                                                                  mRNA.end,
                                                                  len(seq),
                                                                  cds_start_position,
                                                                  cds_end_position)
            cdsRecord = SeqRecord(seq, id=t.replace('transcript:',''), description=desc)
            cds_seq_list.append(cdsRecord)
        else:
            # csv
            seq = str(seq) # for polars write
            it = [g, t, mRNA.chr, mRNA.start, mRNA.end, mRNA.strand, seq, len(seq)]
            if symbol is not None:
                it.insert(1, symbol)
            cds_seq_list.append(it)
    return cds_seq_list

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def split_list(l: list, parts: int) -> list:
    """Takes a list as input, and splits it into "parts" number of sub-lists,
    which are then inserted as elements in the returned meta-list.
    
    The function will try to make the sub-lists as equal in length as 
    possible, so running
    split_list( [1, 2, 3, 4, 5, 6], 4 ) 
    will return 
    [ [1, 2], [3, 4], [5], [6] ]
    
    I will also make sure that the list isn't split into more parts than
    there are elements in the list, so
    split_list( [a, b], 6 ) 
    will return 
    [ [a], [b] ]
    """
    n = min(parts, max(len(l),1))
    k, m = divmod(len(l), n)
    return [l[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]




def get_cds(args):
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
        filename="CDS.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)

    
    if not args.transcript:
        # single cpu
        logging.info("CDS processing at whole-genome scale")
        param_iter = (
            (g, symbol, t, transcript, genome, args.style, args.output_format)
            for g, symbol, t, transcript in iter_transcript_models_with_symbols(
                db, ('CDS', 'exon', 'stop_codon'),
                include_feature_attributes=False,
            )
            if transcript.get('CDS')
        )
        if hasattr(db, 'count_transcripts_with_feature'):
            progress_total = db.count_transcripts_with_feature('CDS')
        else:
            progress_total = sum(
                1 for gene in db.values() for t, transcript in gene.items()
                if t not in {'gene', 'symbol'} and transcript.get('CDS')
            )
        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_CSV_HEADER) + '\n')
        for param in tqdm(param_iter, ncols = 80, total=progress_total, \
                         desc='CDS processing', colour="GREEN"):
            # cds_seq_list.append(sub_cds(param))
            cds_seq_list = sub_cds_faster(param)
            if cds_seq_list == None:
                continue
            stream_parse_output(args, cds_seq_list, out)
        
        # cds_seq_list = [d for de in cds_seq_list if de != None for d in de]
        # parse_output(args, cds_seq_list)
        out.close()
        logging.info("CDS finished at whole-genome scale")
        ## multiprocessing
        """
        param_list = [(g, t, db[g][t], args.genome, args.style, args.output_format) for g in db for t in db[g] if t not in {'gene', 'symbol'} and db[g][t].get('CDS')]
        list_parts = list(np.array_split(param_list, args.process))
        with Pool(processes=args.process) as p:
            cds_seq_list = p.map(multi_mediator, list_parts)
        cds_seq_list = [d[0] for de in cds_seq_list if de != None for d in de]
        """

        """
        # method 2
        queue = Queue()
        for i in range(args.process):
            p = Process(target=multi_mediator2, args=(list_parts[i], queue))
            p.start()
            p.join()
        cds_seq_list = queue.get() 
        """

    else:
        # only one transcript
        logging.info("CDS processing at transcript scale")
        if not t2g.get(args.transcript):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        gene_id = t2g[args.transcript]
        param = (
            gene_id, get_gene_symbol(db, gene_id), args.transcript, db,
            genome, args.style, args.output_format,
        )
        cds_seq_list = sub_cds_faster(param)
        parse_output(args, cds_seq_list)
        logging.info("CDS finished at transcript scale")
