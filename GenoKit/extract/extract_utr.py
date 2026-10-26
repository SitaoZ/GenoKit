# -*- coding: utf-8 -*-
'''
filename: extract_utr.py
date: 2025/11/27 17:07
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import deque
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database, transcript_count,
)
from GenoKit.utils.util import (
    utr3_type, utr5_type, mRNA_type, parse_output, stream_parse_output,
)
from GenoKit.utils.gtf_gff_tool import gff_attr_to_gtf, gtf_attr_to_gff



_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Chrom',
    'Start', 'End', 'Strand',
    'UTR5', 'UTR5 Region', 'UTR5 Count',
    'UTR3', 'UTR3 Region', 'UTR3 Count'
]

def sub_utr(params):
    if len(params) == 8:
        g, t, db, genome, style, output_format, utr3_t, utr5_t = params
        symbol = None
    else:
        g, symbol, t, db, genome, style, output_format, utr3_t, utr5_t = params
    utr_seq_list = deque()
    seq3, seq5 = '', ''
    seq3_count, seq5_count = 0, 0
    seq3_region, seq5_region = [], []
    seq3_gff_lines = deque()
    seq5_gff_lines = deque()
    # utr3
    mRNA = db[g][t]['mRNA']
    if db[g][t].get(utr3_t):
        for u in db[g][t].get(utr3_t):
            seq3_gff_lines.append(u)
            chrom, start, end, strand = u.chr, u.start, u.end, u.strand
            s = str(genome[chrom][start-1 : end])
            seq3 += s
            seq3_count += 1
            seq3_region.append('-'.join(map(str,[u.start,u.end])))
    # utr5
    if db[g][t].get(utr5_t):
        for u in db[g][t].get(utr5_t):
            seq5_gff_lines.append(u)
            chrom, start, end, strand = u.chr, u.start, u.end, u.strand
            s = str(genome[chrom][start-1 : end])
            seq5 += s
            seq5_count += 1
            seq5_region.append('-'.join(map(str,[u.start,u.end])))
    seq3 = Seq(seq3)
    seq5 = Seq(seq5)
    if mRNA.strand == '-':
        seq3 = seq3.reverse_complement()
        seq5 = seq5.reverse_complement()
    if output_format == 'gff':
        for record in seq5_gff_lines:
            # record to gff
            array = list(record)
            if "=" not in array[8]:
                array[8] = gtf_attr_to_gff(array[8])
            gff_line = "\t".join(map(str, array))
            utr_seq_list.append(gff_line)
            # record to gff
            # gff_format = f"{record.chr}\t{record.source}\t{record.type}\t{record.start}\t{record.end}\t{record.score}\t{record.strand}\t{record.phase}\t{record.attribute}"
            # utr_seq_list.append(gff_format)
        for record in seq3_gff_lines:
            array = list(record)
            if "=" not in array[8]:
                array[8] = gtf_attr_to_gff(array[8])
            gff_line = "\t".join(map(str, array))
            utr_seq_list.append(gff_line)
            # gff_format = f"{record.chr}\t{record.source}\t{record.type}\t{record.start}\t{record.end}\t{record.score}\t{record.strand}\t{record.phase}\t{record.attribute}"
            # utr_seq_list.append(gff_format)

    elif output_format == 'gtf':
        for record in seq5_gff_lines:
            # record to gff
            array = list(record)
            if "=" in array[8]:
                array[8] = gff_attr_to_gtf(array[8])
            gff_line = "\t".join(map(str,array))
            utr_seq_list.append(gff_line)
        for record in seq3_gff_lines:
            array = list(record)
            if "=" in array[8]:
                array[8] = gff_attr_to_gtf(array[8])
            gff_line = "\t".join(map(str,array))
            utr_seq_list.append(gff_line)


    elif output_format == 'fasta':
        if len(seq5) != 0:
            utr5_id = ".".join(map(str,[t, '|'.join(seq5_region), '5UTR']))
            desc='strand:%s UTR=5UTR length=%d'%(mRNA.strand, len(seq5))
            utr5Record = SeqRecord(seq5, id=utr5_id.replace('transcript:',''), description=desc)
            utr_seq_list.append(utr5Record)
        if len(seq3) != 0:
            utr3_id = ".".join(map(str,[t, '|'.join(seq3_region), '3UTR']))
            desc='strand:%s UTR=3UTR length=%d'%(mRNA.strand, len(seq3))
            utr3Record = SeqRecord(seq3, id=utr3_id.replace('transcript:',''), description=desc)
            utr_seq_list.append(utr3Record)
    else:
        if len(seq5) != 0 or len(seq3) != 0:
            seq5 = str(seq5) # for polars write
            seq3 = str(seq3)
            it = [g, t.replace('transcript:',''), mRNA.chr, mRNA.start, mRNA.end, mRNA.strand,\
              seq5, '|'.join(seq5_region), seq5_count,\
              seq3, '|'.join(seq3_region), seq3_count]
            header = _CSV_HEADER
            if symbol is not None:
                it.insert(1, symbol)
            else:
                header = [column for column in _CSV_HEADER if column != 'Gene Symbol']
            utr_seq_list.append(dict(zip(header, it)))
    return utr_seq_list
    

def get_utr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        file write to file or stdout
    '''
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Utr.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    # Fasta object
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    utr3_t = utr3_type(args.style)
    utr5_t = utr5_type(args.style)
    if not args.transcript:
        logging.info("UTR processing at whole-genome scale")
        param_iter = (
            (
                g, symbol, t, {g: {t: transcript}}, genome, args.style,
                args.output_format, utr3_t, utr5_t,
            )
            for g, symbol, t, transcript in iter_transcript_models_with_symbols(
                db, (utr3_t, utr5_t),
                include_feature_attributes=args.output_format in {'csv', 'gff', 'gtf'},
            )
        )
        if args.output:
            with open(args.output, 'w') as out:
                if args.output_format == 'csv':
                    out.write(','.join(_CSV_HEADER) + '\n')
                for para in tqdm(
                    param_iter, ncols=80, total=transcript_count(db, t2g),
                    desc='UTR processing:', colour='GREEN',
                ):
                    records = sub_utr(para)
                    if not records:
                        continue
                    if args.output_format == 'csv':
                        for item in records:
                            out.write(','.join(str(item[column]) for column in _CSV_HEADER) + '\n')
                    else:
                        stream_parse_output(args, records, out)
                logging.info("UTR finished at whole-genome scale")
        else:
            # stdout output keeps the historical single-header formatting.
            utr_seq_list = deque()
            for para in tqdm(
                param_iter, ncols=80, total=transcript_count(db, t2g),
                desc='UTR processing:', colour='GREEN',
            ):
                records = sub_utr(para)
                if records:
                    utr_seq_list.extend(records)
            parse_output(args, utr_seq_list)
    else:
        logging.info("UTR processing at transcript scale")
        if not t2g.get(args.transcript):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        # return a specific transcript
        gene_id = t2g[args.transcript]
        params = (
            gene_id, get_gene_symbol(db, gene_id), args.transcript, db,
            genome, args.style, args.output_format, utr3_t, utr5_t,
        )
        utr_seq_list = sub_utr(params)
        # output 
        parse_output(args, utr_seq_list)
        logging.info("UTR finished at transcript scale")



def utr_genbank(args):
    '''
    parameters:
        parse from argparse
    return:
        elements write to a file or stdout
    '''
    genbank = args.genbank
