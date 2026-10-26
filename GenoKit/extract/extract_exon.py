# -*- coding: utf-8 -*-
'''
filename: extract_exon.py
date: 2025/08/12 11:20
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.utils.util import parse_output, stream_parse_output
from GenoKit.database.database import gb_iter, genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database, transcript_count,
)
from GenoKit.utils.gtf_gff_tool import gtf_attr_to_gff, gff_attr_to_gtf

_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID','Chrom','Start','End',
    'Strand', 'Exon', 'Exon Length', 'Exon Index'
]


def _exon_rank(genomic_index, exon_count, strand):
    """Return the 1-based exon number in transcript 5'-to-3' direction."""
    if strand == '-':
        return exon_count - genomic_index + 1
    return genomic_index

def sub_exon(params):
    """
    parameters:
        params:
    return:
        exon_seq_list, a list class
    """
    if len(params) == 6:
        g, t, tdb, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, tdb, genome, style, output_format = params
    exon_seq_list = []
    exon_index = 0
    
    if output_format == 'gff':
        # gff
        for line in tdb['exon']:
            exon_seq_list.append("\t".join(map(str, list(line))))
    else:
        mRNA = tdb['mRNA']
        exon_count = len(tdb['exon'])
        for genomic_index, e in enumerate(tdb['exon'], 1):
            seq = genome[e.chr][e.start-1 : e.end].seq
            seq = Seq(seq)
            if e.strand == '-':
                seq = seq.reverse_complement()
            exon_index = _exon_rank(genomic_index, exon_count, e.strand)
            if output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s exon %s start:%d end:%d length=%d'%(e.strand,
                                                                  exon_index,
                                                                  e.start,
                                                                  e.end,
                                                                  len(seq))
                exonRecord = SeqRecord(seq, id=t, description=desc)
                exon_seq_list.append(exonRecord)
            else:
                # csv
                seq = str(seq) # for polars write
                it = [g, t, e.chr, e.start, e.end, e.strand, seq, len(seq), exon_index]
                header = _CSV_HEADER
                if symbol is not None:
                    it.insert(1, symbol)
                else:
                    header = [column for column in _CSV_HEADER if column != 'Gene Symbol']
                exon_seq_list.append(dict(zip(header, it)))
    return exon_seq_list


def sub_exon_faster(params):
    """
    parameters:
        params:
    return:
        exon_seq_list, a list class
    """
    if len(params) == 6:
        g, t, tdb, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, tdb, genome, style, output_format = params
    exon_seq_list = []
    exon_index = 0

    if output_format == 'gff':
        # gff
        for genomic_index, exon in enumerate(tdb['exon'], 1):
            exon_index = _exon_rank(
                genomic_index, len(tdb['exon']), exon.strand
            )
            array = list(exon)
            if "=" not in array[8]: # gtf to gff
                # array[8] = gtf_attr_to_gff(array[8])
                # Jbrowser2 needs Parent
                array[8] = f"ID={t}_exon{exon_index};Parent={t};GeneID={g};Name={t}_exon{exon_index};"
            # gff to gff no transfer
            gff_line = "\t".join(map(str, array))
            exon_seq_list.append(gff_line)

    elif output_format == 'gtf':
        for genomic_index, exon in enumerate(tdb['exon'], 1):
            exon_index = _exon_rank(
                genomic_index, len(tdb['exon']), exon.strand
            )
            # record to list
            array = list(exon)
            if "=" in array[8]:# gff to gtf
                # array[8] = gff_attr_to_gtf(array[8])
                array[8] = f'exon_id "{t}_exon{exon_index}"; transcript_id "{t}"; gene_id "{g}"; exon_name "{t}_exon{exon_index}";'
            # gtf to gtf no transfer
            gtf_line = "\t".join(map(str, array))
            exon_seq_list.append(gtf_line)

    else:
        mRNA = tdb['mRNA']
        exon_count = len(tdb['exon'])
        for genomic_index, e in enumerate(tdb['exon'], 1):
            seq = genome[e.chr][e.start - 1: e.end].seq
            seq = Seq(seq)
            if e.strand == '-':
                seq = seq.reverse_complement()
            exon_index = _exon_rank(genomic_index, exon_count, e.strand)
            if output_format == 'fasta':
                # fasta (defalut)
                desc = 'strand:%s exon %s start:%d end:%d length=%d' % (e.strand,
                                                                        exon_index,
                                                                        e.start,
                                                                        e.end,
                                                                        len(seq))
                exonRecord = SeqRecord(seq, id=t+f"_{exon_index}", description=desc)
                exon_seq_list.append(exonRecord)
            else:
                # csv
                seq = str(seq)  # for polars write
                it = [g, t, e.chr, e.start, e.end, e.strand, seq, len(seq), exon_index]
                if symbol is not None:
                    it.insert(1, symbol)
                exon_seq_list.append(it)
    return exon_seq_list


def get_exon(args):
    '''
    parameters:
    args.database
    args.output_format
    args.genome
    args.transcript
    args.output
    args.rna_feature
    args.style
    args.print
    '''
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Exon.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    if not args.transcript:
        logging.info("Exon processing at whole-genome scale")
        param_iter = (
            (g, symbol, t, transcript, genome, args.style, args.output_format)
            for g, symbol, t, transcript in iter_transcript_models_with_symbols(
                db, ('exon',),
                include_feature_attributes=args.output_format in {'csv', 'gff', 'gtf'},
            )
        )

        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_CSV_HEADER) + '\n')
        for para in tqdm(param_iter, ncols = 80, total=transcript_count(db, t2g), \
                         desc='Exon processing:', colour="GREEN"):
            exon_seq_list = sub_exon_faster(para)
            if exon_seq_list == None:
                continue
            stream_parse_output(args, exon_seq_list, out)
        
        # exon_seq_list = [d for de in exon_seq_list if de != None for d in de]
        # parse_output(args, exon_seq_list)
        out.close()
        logging.info("Exon finished at whole-genome scale")
    else:
        logging.info("Exon processing at transcript scale")
        t = args.transcript
        if not t2g.get(t):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        g = t2g[t]
        para = (
            g, get_gene_symbol(db, g), t, db[g][t], genome,
            args.style, args.output_format,
        )
        # exon_seq list : 2D list
        exon_seq_list = [sub_exon(para)]
        exon_seq_list = [d for de in exon_seq_list if de != None for d in de]
        # output
        parse_output(args, exon_seq_list)
        logging.info("Exon finished at transcript scale")

def get_exon_gb(args):
    exons = []
    for record in gb_iter(args.genbank):
        index = 1
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                exon_seq = feature.extract(record.seq)
                # feature.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                exon_seq_record = SeqRecord(exon_seq,id='%s exon %d'%(record.id, index),
                                 description='strand %s length %d'%(feature.strand, len(exon_seq)))
                exons.append(exon_seq_record)
                index += 1 
    if args.print:
        SeqIO.write(exons, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(exons, args.output, 'fasta')
 
