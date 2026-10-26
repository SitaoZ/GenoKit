# -*- coding: utf-8 -*-
'''
filename: extract_terminator.py
date: 2024/11/21 15:23
author: Sitao Zhu
mail: zhusitao1990@163.com
'''
import sys
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from collections import deque
from Bio.SeqRecord import SeqRecord
from GenoKit.utils.util import parse_output
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_gene_records_with_symbols, load_annotation_database,
)

# https://pythonhosted.org/gffutils/autodocs/gffutils.Feature.html
# 1-based coordinates

_TERMINATOR_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Chrom', 'Start', 'End',
    'Strand', 'Terminator', 'Length'
]

def sub_terminator(params):
    """
    parameters:
        params: args for multiple processing 
    return:
        deque list 
    """
    terminator_seq_list = deque()
    g, symbol, gdb, genome, terminator_length, utr3_lower_length, output_format = params
    if gdb.strand == "+":
        t_start = gdb.end - utr3_lower_length     # transform 1-based to 0-based
        if t_start < 0:
            t_start = 0
        t_end = gdb.end + int(terminator_length)  # transform 1-based to 0-based
        if t_end > len(genome[gdb.chr]):
            t_end = gdb.end
        terminator = genome[gdb.chr][t_start : t_end].seq # pyfaidx
    else:
        # g.strand == "-"
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].reverse.complement.seq
        t_start = gdb.start - terminator_length - 1  # to 0-based
        if t_start < 0:
            t_start = 0
        t_end = gdb.start + utr3_lower_length - 1 # to 0-based
        if t_end > len(genome[gdb.chr]):
            t_end = gdb.end
        terminator = genome[gdb.chr][t_start : t_end].reverse.complement.seq
    t_start_in_genome = t_start + 1 # back to 1-based 
    t_end_in_genome = t_end


    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d'%(gdb.strand, t_start_in_genome,
                                                 t_end_in_genome, len(terminator))
        terminator = Seq(terminator) # str to Seq
        terminatorRecord = SeqRecord(terminator, id=g, description=desc)
        terminator_seq_list.append(terminatorRecord)

    elif output_format == 'gff':
        # gff
        gff_reacord = [gdb.chr, 'GenoKit', 'terminator', t_start_in_genome, t_end_in_genome,
                           '.', gdb.strand, '.', f'ID={g};']
        terminator_seq_list.append("\t".join(map(str,gff_reacord)))

    elif output_format == 'gtf':
        # gtf
        gtf_reacord = [gdb.chr, 'GenoKit', 'terminator', t_start_in_genome, t_end_in_genome,
                           '.', gdb.strand, '.', f'gene_id "{g}";']
        terminator_seq_list.append("\t".join(map(str,gtf_reacord)))

    else:
        # default csv
        it = [g, symbol, gdb.chr, t_start_in_genome, t_end_in_genome,
              gdb.strand, terminator, len(terminator)]
        terminator_seq_list.append(dict((_TERMINATOR_CSV_HEADER[i], it[i]) for i in range(len(_TERMINATOR_CSV_HEADER))))
    return terminator_seq_list

    

def get_terminator(args):
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
        filename="Terminator.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    terminator_seq_list = deque()
    if not args.gene:
        logging.info("Terminator processing at whole-genome scale")
        # g, db, genome, terminator_length, utr3_lower_length, output_format = params
        param_iter = (
            (g, symbol, gene, genome, args.terminator_length,
             args.utr3_lower_length, args.output_format)
            for g, symbol, gene in iter_gene_records_with_symbols(db)
        )
        for para in tqdm(param_iter, ncols = 80, total=len(db), \
                         desc='Terminator processing:', colour="GREEN"):
            terminator_seq_list.append(sub_terminator(para))
        terminator_seq_list = [d for de in terminator_seq_list if de != None for d in de]
        parse_output(args, terminator_seq_list)
        logging.info("Terminator finished at whole-genome scale")
    else:
        logging.info("Terminator processing at gene scale")
        if not db.get(args.gene):
            print('gene id not exist in database, please check.')
            sys.exit(1)
        g = args.gene
        param = (
            g, get_gene_symbol(db, g), db[g]['gene'], genome,
            args.terminator_length, args.utr3_lower_length, args.output_format,
        )
        terminator_seq_list = sub_terminator(param)
        parse_output(args, terminator_seq_list)
        logging.info("Terminator finished at gene scale")
