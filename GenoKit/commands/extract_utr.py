# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO 
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import genome_dict
from GenoKit.utils.util import utr3_type, utr5_type, mRNA_type, parse_output
from GenoKit.utils.util import gtf_feature_dict, gff_feature_dict



_CSV_HEADER = ['TranscriptID','Chrom', \
               'Start','End','Strand', \
               'UTR5','UTR5_Region','UTR5_Count', \
               'UTR3','UTR3_Region','UTR3_Count']

def sub_utr(params):
    g, t, db, genome, style, output_format, utr3_t, utr5_t = params
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
        for line in seq5_gff_lines:
            utr_seq_list.append(line)
        for line in seq3_gff_lines:
            utr_seq_list.append(line)
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
        it = [t.replace('transcript:',''), mRNA.chr, mRNA.start, mRNA.end, mRNA.strand,\
              seq5, '|'.join(seq5_region), seq5_count,\
              seq3, '|'.join(seq3_region), seq3_count]
        utr_seq_list.append(dict(zip(_CSV_HEADER, it)))
    return utr_seq_list
    

def get_utr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        file write to file or stdout
    '''
    if args.style == 'GFF':
        db, t2g = gff_feature_dict(args.database, args.style)
    else:
        db, t2g = gtf_feature_dict(args.database, args.style)
    genome = genome_dict(args.genome)
    utr3_t = utr3_type(args.style)
    utr5_t = utr5_type(args.style)
    if not args.transcript:
        param_list = [(t2g[t], t, db, genome, args.style, args.output_format, utr3_t, utr5_t) for g in db for t in db[g]]
        utr_seq_list = deque()
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='UTR processing:'):
            utr_seq_list.append(sub_utr(para))
        utr_seq_list = [d for de in utr_seq_list if de != None for d in de]
        parse_output(args, utr_seq_list)
    else:
        # return a specific transcript
        params = (t2g[args.transcript], args.transcript, db, genome, args.style, args.output_format, utr3_t, utr5_t)
        utr_seq_list = sub_utr(params)
        # output 
        parse_output(args, utr_seq_list)



def utr_genbank(args):
    '''
    parameters:
        parse from argparse
    return:
        elements write to a file or stdout
    '''
    genbank = args.genbank
