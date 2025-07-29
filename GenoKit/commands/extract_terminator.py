# -*- coding: utf-8 -*-
import sys
import gffutils
import _pickle as cPickle
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from collections import deque
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.utils.util import parse_output
from GenoKit.database.database import genome_dict

# https://pythonhosted.org/gffutils/autodocs/gffutils.Feature.html
# 1-based coordinates

_CSV_HEADER = ['GeneID','Chrom','Start','End','Strand','Terminator']

def sub_terminator(params):
    """
    parameters:
        params: args for multiple processing 
    return:
        deque list 
    """
    terminator_seq_list = deque()
    g, gdb, genome, terminator_length, utr3_lower_length, output_format = params
    if gdb.strand == "+":
        t_start = gdb.end - utr3_lower_length     # transform 1-based to 0-based
        t_end = gdb.end + int(terminator_length)  # transform 1-based to 0-based
        if t_end > len(genome[gdb.chr]):
            t_end = gdb.end
        terminator = genome[gdb.chr][t_start : t_end].seq # pyfaidx
    else:
        # g.strand == "-"
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].reverse.complement.seq
        t_start = gdb.start - utr3_lower_length - 1 # to 0-based 
        t_end = gdb.start + terminator_length - 1    # to 0-based
        terminator = genome[gdb.chr][t_start : t_end].reverse.complement.seq
    t_start_in_genome = t_start + 1 # back to 1-based 
    t_end_in_genome = t_end
    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d'%(gdb.strand, t_start_in_genome,
                                                 t_end_in_genome, len(terminator))
        terminator = Seq(terminator) # str to Seq
        terminatorRecord = SeqRecord(terminator, id=g, description=desc)
        terminator_seq_list.append(terminatorRecord)
    else:
        # default csv
        it = [g, gdb.chr, t_start_in_genome, \
              t_end_in_genome, gdb.strand, terminator]
        terminator_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return terminator_seq_list

    

def get_terminator(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    terminator_seq_list = deque()
    if not args.gene:
        # g, db, genome, terminator_length, utr3_lower_length, output_format = params
        param_list = [(g, db[g]['gene'], genome, args.terminator_length, args.utr3_lower_length, args.output_format) for g in db]
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='Terminator processing:'):
            terminator_seq_list.append(sub_terminator(para))
        terminator_seq_list = [d for de in terminator_seq_list if de != None for d in de]
        parse_output(args, terminator_seq_list)
    else:
        if not db.get(args.gene):
            sys.exit(1)
            print('gene id not exist in database, please check.')
        g = args.gene
        param = (g, db[g]['gene'], genome, args.terminator_length, args.utr3_lower_length, args.output_format)
        terminator_seq_list = sub_terminator(param)
        parse_output(args, terminator_seq_list)
