# -*- coding: utf-8 -*-
'''
filename: extract_promoter.py
date: 2024/11/21 16:30
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

_PROMOTER_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Chrom', 'Start', 'End',
    'Strand','Promoter', 'Length'
]

_PROMOTER_CSV_HEADER_LINE = ','.join(_PROMOTER_CSV_HEADER) + '\n'
_CSV_WRITE_BUFFER_SIZE = 1024 * 1024


def _promoter_csv_row(gene_id, symbol, gene, genome, promoter_length, utr5_upper_length):
    """Build one promoter CSV row using the established coordinate rules."""
    if gene.strand == '+':
        p_start = gene.start - promoter_length - 1
        if p_start < 0:
            p_start = 0
        p_end = gene.start + utr5_upper_length - 1
        if p_end > len(genome[gene.chr]):
            p_end = len(genome[gene.chr])
        promoter = genome[gene.chr][p_start:p_end].seq
    else:
        p_start = gene.end - utr5_upper_length
        if p_start < 0:
            p_start = 0
        p_end = gene.end + promoter_length
        if p_end > len(genome[gene.chr]):
            p_end = len(genome[gene.chr])
        promoter = genome[gene.chr][p_start:p_end].reverse.complement.seq

    promoter = promoter.upper()
    return (f'{gene_id},{symbol},{gene.chr},{p_start + 1},{p_end},{gene.strand},'
            f'{promoter},{len(promoter)}\n')


def _write_all_promoters_csv(db, genome, args):
    """Write whole-genome promoter CSV without retaining all sequences in RAM."""
    buffered_rows = []
    buffered_size = 0
    with open(args.output, 'w') as out:
        out.write(_PROMOTER_CSV_HEADER_LINE)
        for gene_id, symbol, gene in tqdm(iter_gene_records_with_symbols(db), ncols=80, total=len(db),
                                       desc='Promoter processing:', colour='GREEN'):
            row = _promoter_csv_row(
                gene_id, symbol, gene, genome, args.promoter_length,
                args.utr5_upper_length)
            buffered_rows.append(row)
            buffered_size += len(row)
            if buffered_size >= _CSV_WRITE_BUFFER_SIZE:
                out.write(''.join(buffered_rows))
                buffered_rows.clear()
                buffered_size = 0
        if buffered_rows:
            out.write(''.join(buffered_rows))

def sub_promoter(params):
    """
    parameters:
        params: args for multiple processing 
    return:
        deque list 
    """
    promoter_seq_list = deque()
    g, symbol, gdb, genome, promoter_length, utr5_upper_length, output_format = params
    if gdb.strand == "+":
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].seq
        p_start = gdb.start - promoter_length - 1 # 往前数,g.end不需要, g.start需要-1, transform 1-based to 0-based
        if p_start < 0 :
            p_start = 0
        p_end = gdb.start + utr5_upper_length - 1 # include first base of TSS, transform 1-based to 0-based
        if p_end > len(genome[gdb.chr]):
            p_end = len(genome[gdb.chr])
        promoter = genome[gdb.chr][p_start : p_end].seq # pyfaidx
    else:
        # g.strand == "-"
        # gene_seq = genome[gdb.chr][gdb.start : gdb.end].reverse.complement.seq
        p_start = gdb.end - utr5_upper_length
        if p_start < 0 :
            p_start = 0
        p_end = gdb.end + promoter_length
        if p_end > len(genome[gdb.chr]):
            p_end = len(genome[gdb.chr])
        promoter = genome[gdb.chr][p_start : p_end].reverse.complement.seq

    promoter = promoter.upper()
    p_start_in_genome = p_start + 1 # back to 1-based 
    p_end_in_genome = p_end

    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d'%(gdb.strand, p_start_in_genome,
                                                 p_end_in_genome, len(promoter))
        promoter = Seq(promoter) # str to Seq
        promoterRecord = SeqRecord(promoter, id=g, description=desc)
        promoter_seq_list.append(promoterRecord)
    elif output_format == 'gff':
        # gff
        gff_reacord = [gdb.chr, 'GenoKit', 'promoter', p_start_in_genome, p_end_in_genome,
                       '.', gdb.strand, '.', f'ID={g};']
        promoter_seq_list.append("\t".join(map(str, gff_reacord)))
    elif output_format == 'gtf':
        # gtf
        gtf_reacord = [gdb.chr, 'GenoKit', 'promoter', p_start_in_genome, p_end_in_genome,
                       '.', gdb.strand, '.', f'gene_id "{g}";']
        promoter_seq_list.append("\t".join(map(str, gtf_reacord)))
    else:
        # default csv
        it = [g, symbol, gdb.chr, p_start_in_genome, \
              p_end_in_genome, gdb.strand, promoter, len(promoter)]
        promoter_seq_list.append(dict((_PROMOTER_CSV_HEADER[i], it[i]) for i in range(len(_PROMOTER_CSV_HEADER))))
    return promoter_seq_list

    
def get_promoter(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Promoter.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    promoter_seq_list = deque()
    if not args.gene:
        logging.info("Promoter processing at whole-genome scale")
        if args.output and args.output_format == 'csv':
            _write_all_promoters_csv(db, genome, args)
            return

        # g, db, genome, promoter_length, utr5_upper_length, output_format = params
        param_iter = (
            (g, symbol, gene, genome, args.promoter_length,
             args.utr5_upper_length, args.output_format)
            for g, symbol, gene in iter_gene_records_with_symbols(db)
        )
        for para in tqdm(param_iter, ncols = 80, total=len(db), \
                         desc='Promoter processing:', colour="GREEN"):
            promoter_seq_list.append(sub_promoter(para))
        promoter_seq_list = [d for de in promoter_seq_list if de != None for d in de]
        parse_output(args, promoter_seq_list)
        logging.info("Promoter finished at whole-genome scale")
    else:
        logging.info("Promoter processing at gene scale")
        if not db.get(args.gene):
            print('gene id not exist in database, please check.')
            sys.exit(1)
        g = args.gene
        param = (
            g, get_gene_symbol(db, g), db[g]['gene'], genome,
            args.promoter_length, args.utr5_upper_length, args.output_format,
        )
        promoter_seq_list = sub_promoter(param)
        parse_output(args, promoter_seq_list)
        logging.info("Promoter finished at gene scale")
