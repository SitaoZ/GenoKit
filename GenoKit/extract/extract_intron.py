#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: extract_intron.py
date: 2025/12/18 13:34
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
import itertools
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database, transcript_count,
)
from GenoKit.utils.util import parse_output, stream_parse_output



_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Intron', 'Chrom',
    'Start Genome', 'End Genome', 'Strand', 'Sequence'
]


def sub_intron(params):
    '''
    parameters:
        params: argparse
    return:
        intron_seq_list, a deque
    '''

    # g, t, db[g][t], genome, args.style, args.output_format
    if len(params) == 6:
        g, t, tdb, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, tdb, genome, style, output_format = params
    mRNA = tdb['mRNA']
    # intron_seq_list = deque()
    intron_seq_list = []
    exons = list(itertools.chain(*[[e.start,e.end] for e in tdb['exon']])) # exon position include start and end
    if len(exons)/2 == 1:
        # only one exon
        return None
    else:
        left = exons[1:-1]
        intron_index = 0
        for i in range(0, len(left), 2):
            s = left[i]   # intron start
            e = left[i+1] # intron end, shoukd be minus 1 
            if s+1 != e:
                # exon 与exon 不相连，就是intron的位置
                # 注意索引的点，
                intron_index += 1
                intron_seq= genome[mRNA.chr][s:e-1]
                if mRNA.strand == '-':
                    intron_seq = intron_seq.reverse
                    intron_seq = intron_seq.complement
                    intron_rank = int(len(exons)/2 - intron_index)
                else:
                    intron_rank = intron_index
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    intron_gff = [mRNA.chr, 'GenoKit', 'intron', s+1, 
                                  e-1, '.', mRNA.strand, '.', f'ID={t};Name=intron{intron_rank};']
                    intron_gff_line = list(map(str, intron_gff))
                    intron_seq_list.append("\t".join(intron_gff_line))
                    continue
                elif output_format == 'fasta':
                    intron_seq= Seq(str(intron_seq))
                    # 1 based position
                    intronRecord = SeqRecord(intron_seq, id=t, 
                                             description='strand %s intron%d start %d end %d length=%d exons=%d'%(
                                             mRNA.strand, intron_rank, s+1, e-1, len(intron_seq), len(exons)/2)
                                            )
                    intron_seq_list.append(intronRecord)
                else:
                    intron_seq= str(intron_seq) # for polars write
                    it = [g, t, f'intron{intron_rank}', mRNA.chr, s+1, e-1, mRNA.strand, intron_seq]
                    header = _CSV_HEADER
                    if symbol is not None:
                        it.insert(1, symbol)
                    else:
                        header = [column for column in _CSV_HEADER if column != 'Gene Symbol']
                    intron_seq_list.append(dict(zip(header, it)))
                    # intron_seq_list.append(it)
        return intron_seq_list


def sub_intron_faster(params):
    '''
    parameters:
        params:
    return:
        intron_seq_list, a deque
    '''

    # g, t, db[g][t], genome, args.style, args.output_format
    if len(params) == 6:
        g, t, tdb, genome, style, output_format = params
        symbol = None
    else:
        g, symbol, t, tdb, genome, style, output_format = params
    mRNA = tdb['mRNA']
    # intron_seq_list = deque()
    intron_seq_list = []
    exons = list(itertools.chain(*[[e.start,e.end] for e in tdb['exon']])) # exon position include start and end
    if len(exons)/2 == 1:
        # only one exon
        return None
    else:
        left = exons[1:-1]
        intron_index = 0
        for i in range(0, len(left), 2):
            s = left[i]   # intron start
            e = left[i+1] # intron end, shoukd be minus 1
            if s+1 != e:
                # exon 与exon 不相连，就是intron的位置
                # 注意索引的点，
                intron_index += 1
                intron_seq= genome[mRNA.chr][s:e-1]
                if mRNA.strand == '-':
                    intron_seq = intron_seq.reverse
                    intron_seq = intron_seq.complement
                    intron_rank = int(len(exons)/2 - intron_index)
                else:
                    intron_rank = intron_index
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    intron_gff = [mRNA.chr, 'GenoKit', 'intron', s + 1, e - 1, '.',
                                  mRNA.strand, '.',
                                  f'ID={t}_intron{intron_rank};TranscriptID={t};GeneID={g};Name={t}_intron{intron_rank};']
                    intron_gff_line = list(map(str, intron_gff))
                    intron_seq_list.append("\t".join(intron_gff_line))
                    continue
                elif output_format == 'gtf':
                    intron_gtf = [mRNA.chr, 'GenoKit', 'intron', s + 1, e - 1, '.',
                                  mRNA.strand, '.',
                                  f'intron_id "{t}_intron{intron_rank}"; transcript_id "{t}"; gene_id "{g}"; intron_name "{t}_intron{intron_rank}";']
                    intron_gtf_line = list(map(str, intron_gtf))
                    intron_seq_list.append("\t".join(intron_gtf_line))
                elif output_format == 'fasta':
                    intron_seq= Seq(str(intron_seq))
                    # 1 based position
                    intronRecord = SeqRecord(intron_seq, id=t,
                                             description='strand %s intron%d start %d end %d length=%d exons=%d'%(
                                             mRNA.strand, intron_rank, s+1, e-1, len(intron_seq), len(exons)/2)
                                            )
                    intron_seq_list.append(intronRecord)
                else:
                    intron_seq= str(intron_seq) # for polars write
                    it = [g, t, f'intron{intron_rank}', mRNA.chr, s+1, e-1, mRNA.strand, intron_seq]
                    if symbol is not None:
                        it.insert(1, symbol)
                    intron_seq_list.append(it)
        return intron_seq_list

def get_intron(args):
    '''
    intron file size is big
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Intron.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)
    if not args.transcript:
        logging.info("Intron processing at whole-genome scale")
        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_CSV_HEADER)+'\n')
        param_iter = (
            (g, symbol, t, transcript, genome, args.style, args.output_format)
            for g, symbol, t, transcript in iter_transcript_models_with_symbols(
                db, ('exon',), include_feature_attributes=False
            )
        )
        for para in tqdm(param_iter, ncols = 80, total=transcript_count(db, t2g), \
                         desc='Intron processing:', colour="GREEN"):
            intron_seq_list = sub_intron_faster(para)
            if intron_seq_list == None:
                continue
            stream_parse_output(args, intron_seq_list, out)
        out.close()
        logging.info("Intron finished at whole-genome scale")
    else:
        logging.info("Intron processing at transcript scale")
        t = args.transcript
        if not t2g.get(t):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        g = t2g[t]
        para = (
            g, get_gene_symbol(db, g), t, db[g][t], genome,
            args.style, args.output_format,
        )
        # intron_seq_list : 2D list
        intron_seq_list = [sub_intron(para)]
        if intron_seq_list[0] == None :
            print('The transcript does not have intron, because of one exon.')
            sys.exit(1)
        intron_seq_list = [d for de in intron_seq_list if de != None for d in de]
        parse_output(args, intron_seq_list)
        logging.info("Intron finished at transcript scale")
