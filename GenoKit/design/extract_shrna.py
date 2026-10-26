#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: extract_shrna.py
date: 2025/12/24 20:34
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

# !/usr/bin/env python3
"""
shRNA设计程序 - 用于RNA干扰实验的高效shRNA设计
shRNA 靶向成熟的转录本
特性：
1. 多参数评分系统（热力学、GC含量、连续碱基等）
2. 特异性评估（避免脱靶）
3. 二级结构预测
4. 批量设计和高通量筛选
5. 结果自动保存为CSV文件
"""

import sys
import re
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from typing import Dict, Tuple
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.utils.util import parse_output
from GenoKit.utils.util import stream_parse_output
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import get_gene_symbol, load_annotation_database
from GenoKit.database.kmer_database import open_kmer_database
from GenoKit.design.design_shrna import shRNADesigner


_shRNA_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'shRNA', 'GC', 'Strand', 'Score',
    'Targets', 'Kmer Freq', 'Loop', 'Sense Oligo', 'Antisense Oligo'
]

def _parse_target_sequences(params: Tuple, region: str) -> Dict[str, Dict[str, str]]:
    """Build selected sequences in mature 5'-to-3' orientation."""
    gene_isoform = defaultdict(dict)
    g, transcript_ids, db, genome, _style, _output_format = params
    for transcript_id in transcript_ids:
        transcript = db[g][transcript_id]
        features = transcript['exon']
        if region == 'cds':
            features = transcript.get('CDS') or features
        sequence = ''.join(
            str(genome[feature.chr][feature.start - 1:feature.end])
            for feature in features
        )
        if transcript['mRNA'].strand == '-':
            sequence = str(Seq(sequence).reverse_complement())
        gene_isoform[g][transcript_id] = sequence.upper()
    return gene_isoform


def parse_transcript(params: Tuple) -> Dict[str, Dict[str, str]]:
    """Build spliced mature transcripts, including UTRs."""
    return _parse_target_sequences(params, 'transcript')


def parse_cds(params: Tuple) -> Dict[str, Dict[str, str]]:
    """Build spliced CDS; fall back to exons when CDS is not annotated."""
    return _parse_target_sequences(params, 'cds')


def get_shrna(args):
    '''

    :param args:
    :return:
    '''

    # logging
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="shRNA.log"
    )

    target_mode = getattr(args, 'target_mode', 'all')
    region = getattr(args, 'region', 'cds')
    if region not in {'cds', 'transcript'}:
        raise ValueError("region must be 'cds' or 'transcript'")
    sequence_parser = parse_transcript if region == 'transcript' else parse_cds
    db, t2g = load_annotation_database(args.database)
    genome = genome_dict(args.genome)

    # shRNAs bind mature mRNA, so only transcriptome-labelled catalogues are
    # accepted when source metadata is available.
    kmer_bin_database = open_kmer_database(
        args.kmer_bin,
        k=args.length,
        expected_source='transcriptome',
        jellyfish=getattr(args, 'jellyfish', None),
    )
    logging.info(f"Perform shRNA design at the whole-genome scale")
    if not args.gene:
        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_shRNA_CSV_HEADER) + '\n')
        for g in tqdm(db, ncols=80, total=len(db),
                      desc = 'shRNA processing:', colour="GREEN"):
            gene = db[g]['gene']
            symbol = get_gene_symbol(db, g)
            if isinstance(gene, dict):
                print(gene, g)
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            # param setting
            ts = [i for i in db[g] if i not in {'gene', 'symbol'}]
            params = (g, ts, db, genome, args.style, args.output_format)
            gene_isoform = sequence_parser(params)
            # print(gene_isoform)
            # gene_isoform = parse_cds_exon(params)
            test_transcripts = gene_isoform[g]
            if len(test_transcripts) == 0:
                continue
            designer = shRNADesigner(
                g, test_transcripts, strand,
                kmer_bin=kmer_bin_database,
                shrna_length=args.length,
                loop_type=args.loop,
            )
            shrna_list = designer.design(target_mode=target_mode)

            # fasta (defalut)
            if not shrna_list:
                # no shRNA
                continue
            shrna_seq_list = []
            for index,shrna in enumerate(shrna_list):
                # shrna is a dict
                sequence = shrna['sequence']
                strand = shrna['strand'] # not same to gene.strand
                gc = shrna['gc_content']
                score = shrna['score']
                targets = shrna['targets']
                loop = shrna['loop']
                kmer_freq = shrna['kmer_freq']
                positions = shrna['positions']
                sense_oligo = shrna['sense_oligo']
                antisense_oligo = shrna['antisense_oligo']


                if args.output_format == 'fasta':
                    desc='strand:%s transcript:%s Loop:%s start:%d end:%d length:%d score:%d kmer_freq:%d'%(gene.strand,
                                                                                                            "|".join(t for t in ts if t in positions),
                                                                                                            loop,
                                                                                                            gene.start,
                                                                                                            gene.end,
                                                                                                            len(re.sub(r'[^ATCGatcg]', '', sequence)),
                                                                                                            score,
                                                                                                           kmer_freq)
                    sgRecord = SeqRecord(Seq(sequence), id=g+f"_{index+1}", description=desc)
                    shrna_seq_list.append(sgRecord)
                else:
                    it = [g, symbol, "|".join(t for t in ts if t in positions), sequence, gc, strand, score, targets, kmer_freq, loop, sense_oligo, antisense_oligo]
                    shrna_seq_list.append(it)
            stream_parse_output(args, shrna_seq_list, out)
        # parse_output(args, shrna_seq_list)
        out.close()
        logging.info(f"Finished shRNA design at the whole-genome scale")
    else:
        logging.info(f"Perform shRNA design ath the gene level")
        shrna_seq_list = []
        if args.gene not in db:
            print(f'Error: {args.gene} is not in the annotation database.')
            return None
        gene = db[args.gene]['gene']
        g = args.gene
        symbol = get_gene_symbol(db, g)
        chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
        ts = [i for i in db[args.gene] if i not in {'gene', 'symbol'}]
        params = (args.gene, ts, db, genome, args.style, args.output_format)
        gene_isoform = sequence_parser(params)
        test_transcripts = gene_isoform[args.gene]
        # Apply target-mode to the selected region of each transcript.
        designer = shRNADesigner(
            g, test_transcripts, strand, kmer_bin=kmer_bin_database,
            shrna_length=args.length, loop_type=args.loop,
        )
        shrna_list = designer.design(target_mode=target_mode)
        for shrna in shrna_list:
            sequence = shrna['sequence']
            strand = shrna['strand']  # not same to gene.strand
            gc = shrna['gc_content']
            score = shrna['score']
            targets = shrna['targets']
            loop = shrna['loop']
            kmer_freq = shrna['kmer_freq']
            positions = shrna['positions']
            sense_oligo = shrna['sense_oligo']
            antisense_oligo = shrna['antisense_oligo']

            if args.output_format == 'fasta':
                desc = 'strand:%s transcript:%s PAM:%s start:%d end:%d length:%d score:%d kmer_freq:%d' % (gene.strand,
                                                                                                              "|".join(t for t in ts if t in positions),
                                                                                                              loop,
                                                                                                              gene.start,
                                                                                                              gene.end,
                                                                                                              len(re.sub(r'[^ATCGatcg]', '', sequence)),
                                                                                                              score,
                                                                                                           kmer_freq)
                sgRecord = SeqRecord(Seq(sequence), id=args.gene, description=desc)
                shrna_seq_list.append(sgRecord)
            else:
                it = [args.gene, symbol, "|".join(t for t in ts if t in positions), sequence, gc, strand, score, targets, kmer_freq, loop, sense_oligo, antisense_oligo]
                shrna_seq_list.append(dict((_shRNA_CSV_HEADER[i], it[i]) for i in range(len(_shRNA_CSV_HEADER))))
        parse_output(args, shrna_seq_list)
        logging.info(f"Finished shRNA design at the gene scale")


if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": "ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_002": "ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_003": "ATGAGATCCTGGCCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA"}
    designer = shRNADesigner("TP53", test_transcripts)
    print(designer)
    # 寻找通用sgRNA
    universal = designer.design(require_all=True)
    print(f"找到 {len(universal)} 条通用sgRNA:")
    for sh in universal[:3]:
        print(f"Seq: {sh['sequence']} ({sh['strand']}链)")
        print(f"GC: {sh['gc_content']}, 评分: {sh['score']}")
        print(f"靶向转录本: {len(sh['positions'])}个\n")

    # 寻找任意sgRNA
    any_sh = designer.design(require_all=False)
    print(f"\n总候选sgRNA: {len(any_sh)} 条")
