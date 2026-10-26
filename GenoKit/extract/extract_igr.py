# -*- coding: utf-8 -*-
'''
filename: extract_igr.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
import sys
import logging
import itertools
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.utils.util import parse_output, stream_parse_output
from GenoKit.database.database import gb_iter, genome_dict
from GenoKit.database.annotation_db import (
    iter_gene_records_with_symbols, load_annotation_database,
)


_CSV_HEADER = [
    'IGR ID', 'Gene Symbol', 'Chrom', 'Strand',
    'Start Genome', 'End Genome', 'Sequence'
]

def sub_IGR(params):
    """
    parameters:
        params:
    return:
        igr_seq_list, a list class
    """
    chrom, database, genome, style, output_format, genes, gene_symbols = params
    if style.lower() == 'gff':
        g_pat = re.compile(r'ID=(\S+?);')
    else:
        g_pat = re.compile(r'gene_id "(\S+?)";')
    igr_seq_list = []
    # ['start1', 'end1', 'start2', 'end2', 'start3', 'end3', 'start4', 'end4', 'start5', 'end5']
    gene_positions = list(itertools.chain(*[[g.start,g.end] for g in genes]))
    gene_ids = []
    for g in genes:
        match = g_pat.search(g.attribute)
        gene_ids.append(match.group(1)) 
    strands = [g.strand for g in genes]

    if len(gene_positions)%2 == 1:
        # only one gene
        return None
    else:
        left = gene_positions[1:-1] # index 1
        for i in range(0, len(left), 2):
            # step=2 allows skipping the gene region
            s = left[i]   # igr start
            e = left[i+1] # igr end
            igr_id = gene_ids[int(i/2)]+"_"+gene_ids[int(i/2)+1]
            igr_symbol = gene_symbols[int(i/2)]+"|"+gene_symbols[int(i/2)+1]
            igr_strand = strands[int(i/2)]+"|"+strands[int(i/2)+1]
            if s+1 != e:
                # gene 与gene 不相连，就是igr的位置
                # 注意索引的点，
                if s > e-1:
                    # gene overlap
                    continue 
                igr_seq = genome[chrom][s:e-1] # e should be minus 1
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    igr_gff = [chrom, 'GenoKit', 'igr', s+1, e-1, '.', igr_strand, '.', f'ID={igr_id};']
                    igr_gff_line = list(map(str, igr_gff))
                    igr_seq_list.append("\t".join(igr_gff_line))
                    continue
                elif output_format == 'fasta':
                    igr_seq= Seq(str(igr_seq))
                    # 1 based position
                    igrRecord = SeqRecord(igr_seq, id=f'{igr_id}',
                                             description='strand %s start %d end %d length=%d'%(
                                             igr_strand, s+1, e-1, len(igr_seq))
                                            )
                    igr_seq_list.append(igrRecord)
                else:
                    igr_seq= str(igr_seq)
                    it = [igr_id, igr_symbol, chrom, igr_strand, s+1, e-1, igr_seq]
                    igr_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        return igr_seq_list



def sub_IGR_faster(params):
    """
    parameters:
        params:
    return:
        igr_seq_list, a list class
    """
    chrom, database, genome, style, output_format, genes, gene_symbols = params
    if style.lower() == 'gff':
        g_pat = re.compile(r'ID=(\S+?);')
    else:
        g_pat = re.compile(r'gene_id "(\S+?)";')
    igr_seq_list = []
    # ['start1', 'end1', 'start2', 'end2', 'start3', 'end3', 'start4', 'end4', 'start5', 'end5']
    gene_positions = list(itertools.chain(*[[g.start,g.end] for g in genes]))
    gene_ids = []
    for g in genes:
        match = g_pat.search(g.attribute)
        gene_ids.append(match.group(1))
    strands = [g.strand for g in genes]

    if len(gene_positions)%2 == 1:
        # only one gene
        return None
    else:
        left = gene_positions[1:-1] # index 1
        for i in range(0, len(left), 2):
            # step=2 allows skipping the gene region
            s = left[i]   # igr start
            e = left[i+1] # igr end
            igr_id = gene_ids[int(i/2)]+"_"+gene_ids[int(i/2)+1]
            igr_symbol = gene_symbols[int(i/2)]+"|"+gene_symbols[int(i/2)+1]
            igr_strand = strands[int(i/2)]+"|"+strands[int(i/2)+1]
            if s+1 != e:
                # gene 与gene 不相连，就是igr的位置
                # 注意索引的点，
                if s > e-1:
                    # gene overlap
                    continue
                igr_seq = genome[chrom][s:e-1] # e should be minus 1
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    igr_gff = [chrom, 'GenoKit', 'igr', s+1, e-1, '.', igr_strand, '.', f'ID={igr_id};']
                    igr_gff_line = list(map(str, igr_gff))
                    igr_seq_list.append("\t".join(igr_gff_line))
                    continue
                elif output_format == 'gtf':
                    igr_gtf = [chrom, 'GenoKit', 'igr', s + 1, e - 1, '.', igr_strand, '.', f'gene_id "{igr_id}";']
                    igr_gtf_line = list(map(str, igr_gtf))
                    igr_seq_list.append("\t".join(igr_gtf_line))
                elif output_format == 'fasta':
                    igr_seq= Seq(str(igr_seq))
                    # 1 based position
                    igrRecord = SeqRecord(igr_seq, id=f'{igr_id}',
                                             description='strand %s start %d end %d length=%d'%(
                                             igr_strand, s+1, e-1, len(igr_seq))
                                            )
                    igr_seq_list.append(igrRecord)
                else:
                    igr_seq= str(igr_seq)
                    it = [igr_id, igr_symbol, chrom, igr_strand, s+1, e-1, igr_seq]
                    igr_seq_list.append(it)
        return igr_seq_list


def get_igr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # GTF noncoding genes are annotated transcript, while GFF noncoding genes are annotated gene
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="IGR.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    chrom_split = dict()
    for _gene_id, symbol, gene in iter_gene_records_with_symbols(db):
        chrom = gene.chr
        if chrom_split.get(chrom):
            chrom_split[chrom].append((gene, symbol))
        else:
            chrom_split[chrom] = [(gene, symbol)]
    # sort by coordinate
    for chrom in chrom_split:
        chrom_split[chrom] = sorted(
            chrom_split[chrom], key=lambda item: item[0].start,
        )

    param_list = [
        (
            chrom, args.database, genome, args.style, args.output_format,
            [item[0] for item in chrom_split[chrom]],
            [item[1] for item in chrom_split[chrom]],
        )
        for chrom in chrom_split
    ]
    out = open(args.output, 'w')
    if args.output_format == 'csv':
        out.write(",".join(_CSV_HEADER) + '\n')
    logging.info("IGR processing at whole-genome scale")
    for param in tqdm(param_list, total=len(param_list), ncols=80, \
                      desc='IGR processing:', colour="GREEN"):
        igr_seq_list = sub_IGR_faster(param)
        if igr_seq_list == None:
            continue
        stream_parse_output(args, igr_seq_list, out)
    out.close()
    logging.info("IGR finished at whole-genome scale")

    # igr_seq_list = [d for de in igr_seq_list if de != None for d in de]
    # parse_output(args, igr_seq_list)

def get_IGR_gb(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    gene = []
    for record in gb_iter(args.genbank):
        for feature in record.features:
            if feature.type == 'gene': # CDS promoter UTR 
                gene_seq = ''
                for part in feature.location.parts:
                    gene_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                gene_seq_record = SeqRecord(gene_seq, id='gene:%s'%(gene_id), 
                                  description='strand %s length %d'%(feature.strand, len(gene_seq)))
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
