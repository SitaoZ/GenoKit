# -*- coding: utf-8 -*-
import re
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
import _pickle as cPickle
from tqdm import tqdm
from Bio.Seq import Seq
from collections import deque
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.utils.util import parse_output
from GenoKit.database.database import genome_dict

# https://pythonhosted.org/gffutils/autodocs/gffutils.Feature.html
# 1-based coordinates

_CSV_HEADER = ['GeneID','Chrom','Start','End','Strand','Promoter']

def sub_promoter(params):
    """
    parameters:
        params: args for multiple processing 
    return:
        deque list 
    """
    promoter_seq_list = deque()
    g, gdb, genome, promoter_length, utr5_upper_length, output_format = params
    if gdb.strand == "+":
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].seq
        p_start = gdb.start - promoter_length - 1 # 往前数,g.end不需要, g.start需要-1, transform 1-based to 0-based
        if p_start < 0 :
            p_start = gdb.start - 1
        p_end = gdb.start + utr5_upper_length - 1 # include first base of TSS, transform 1-based to 0-based
        promoter = genome[gdb.chr][p_start : p_end].seq # pyfaidx
    else:
        # g.strand == "-"
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].reverse.complement.seq
        p_start = gdb.end - utr5_upper_length 
        p_end = gdb.end + promoter_length
        promoter = genome[gdb.chr][p_start : p_end].reverse.complement.seq
    p_start_in_genome = p_start + 1 # back to 1-based 
    p_end_in_genome = p_end
    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d'%(gdb.strand, p_start_in_genome,
                                                 p_end_in_genome, len(promoter))
        promoter = Seq(promoter) # str to Seq
        promoterRecord = SeqRecord(promoter, id=g, description=desc)
        promoter_seq_list.append(promoterRecord)
    else:
        # default csv
        it = [g, gdb.chr, p_start_in_genome, \
              p_end_in_genome, gdb.strand, promoter]
        promoter_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return promoter_seq_list

    
def get_promoter(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    promoter_seq_list = deque()
    if not args.gene:
        # g, db, genome, promoter_length, utr5_upper_length, output_format = params
        param_list = [(g, db[g]['gene'], genome, args.promoter_length, args.utr5_upper_length, args.output_format) for g in db]
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='Promoter processing:'):
            promoter_seq_list.append(sub_promoter(para))
        promoter_seq_list = [d for de in promoter_seq_list if de != None for d in de]
        parse_output(args, promoter_seq_list)
    else:
        if not db.get(args.gene):
            sys.exit(1)
            print('gene id not exist in database, please check.')
        g = args.gene
        param = (g, db[g]['gene'], genome, args.promoter_length, args.utr5_upper_length, args.output_format)
        promoter_seq_list = sub_promoter(param)
        parse_output(args, promoter_seq_list)
