# -*- coding: utf-8 -*-
'''
filename: extract_sgrna.py
date: 2025/12/5 11:05
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
import logging
from tqdm import tqdm
from Bio.Seq import Seq
from tabulate import tabulate
from typing import Dict, Tuple
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import get_gene_symbol, load_annotation_database
from GenoKit.database.kmer_database import open_kmer_database
from GenoKit.utils.util import parse_output
from GenoKit.design.design_sgrna import CDSsgRNADesigner

_sgRNA_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'sgRNA', 'GC',
    'Strand', 'Score', 'Targets', 'PAM', 'Kmer Freq'
]
# Targets: 靶向基因内几条转录本

# 设计原则
# 1.sgRNA靶向序列必须紧邻一个特定的原间隔区邻近基序（Protospacer Adjacent Motif, PAM）序列。
# 对于最常用的SpCas9核酸酶，其PAM序列特征为NGG（其中N可以是任意核苷酸，位于靶向序列的3'端）。
# 设计时，必须确保目标位点存在兼容的PAM序列。

# 2.靶向序列（Spacer sequence）：最佳长度通常为20个核苷酸（nt）。
# 序列过短可能会降低编辑效率和特异性；序列过长则会显著降低效率。

# 3.避免脱靶
# 序列的特异性
# 在编码区（CDS）选择靶点: 为了实现基因敲除（Knockout, KO），
# 靶点通常选择在基因编码区的 前 50%，尤其是关键的 功能结构域 附近，
# 以确保产生的截短蛋白（如果翻译发生）是无功能的。
# 避免 ATG 所在外显子: 避免在起始密码子 ATG 所在的外显子或其上游外显子进行敲除，
# 因为这可能会导致基因表达完全丧失，难以获得稳定的敲除细胞株。
# 创建移码突变: 理想情况下，sgRNA 靶点导致的 DNA 双链断裂（DSB）修复（通过 NHEJ 途径）
# 应在编码区引入非三的倍数碱基的插入或缺失，从而产生 移码突变（Frameshift Mutation），
# 导致下游提前终止，高效实现功能性敲除.
#  4.活性，效率
#  GC含量：40%-80%，有助于维持sgRNA:DNA双链的稳定性
#  避免特定的基序：避免sgRNA3'端出现连续4个以上的T
# 5. 目标区域选择
#  基因敲除, 靶向基因的外显子，且编码蛋白的区域
#  基因激活、抑制，靶向转录起始位点

def pam_table():
    # Size is Cas protein length
    # Location is PAM location opposite to cutting site
    PAM = [
        ["Type",     "Size",  "PAM",      "Size of sgRNA", "Location",       "Target"],
        ["SpCas9",   "1368",  "NGG",      "20",            "downstream(3')", "dsDNA"],
        ["SaCas9",   "1053",  "NNGRRT",   "21",            "downstream(3')", "dsDNA"],
        ["CjCas9",   "984",   "NNNNRYAC", "22",            "downstream(3')", "DsDNA"],
        ["St1Cas9",  "1121",  "NNAGAAW",  "20",            "downstream(3')", "dsDNA"],
        ["AsCas12a", "1307",  "TTTV",     "24",            "upstream(5')",   "ssDNA,dsDNA"],
        ["LbCas12a", "1228",  "TTTV",     "24",            "upstream(5')",   "ssDNA,dsDNA"],
        ["Cas12b",   "/",    "TTN",      "/",             "upstream(5')",   "ssDNA,dsDNA"],
        ["Cas14b",   "400~700",    "TTTA",     "/",             "upstream(5')",   "ssDNA"],
    ]
    print(tabulate(PAM, headers="firstrow", tablefmt="grid"))

def parse_transcript(params:Tuple) -> Dict[str, Dict]:
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    for t in ts:
        seq = ''
        mRNA = db[g][t]['mRNA']
        for e in db[g][t]['exon']:
            chrom, start, end, strand = e.chr, e.start, e.end, e.strand
            s = str(genome[chrom][start-1:end])
            seq += s
        seq = Seq(seq)
        if mRNA.strand == '-':
            seq = seq.reverse_complement()
        seq = str(seq).upper()
        gene_isoform[g][t] = seq
    return gene_isoform

def parse_cds(params:Tuple) -> Dict[str, Dict]:
    """
    sgRNA 需要靶向DNA，编码CDS的外显子效果最佳，intron效果不好
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    for t in ts:
        seq = ''
        mRNA = db[g][t]['mRNA']
        if db[g][t].get("CDS"):
            for c in db[g][t]['CDS']:
                chrom, start, end, strand = c.chr, c.start, c.end, c.strand
                s = str(genome[chrom][start-1:end])
                seq += s
            seq = Seq(seq)
            if mRNA.strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()
            gene_isoform[g][t] = seq
        else:
            for e in db[g][t]['exon']:
                chrom, start, end, strand = e.chr, e.start, e.end, e.strand
                s = str(genome[chrom][start - 1:end])
                seq += s
            seq = Seq(seq)
            if mRNA.strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()
            gene_isoform[g][t] = seq
    return gene_isoform

def parse_cds_exon(params:Tuple) -> Dict[str, Dict]:
    """
    sgRNA 需要靶向DNA，编码CDS的外显子效果最佳，intron效果不好
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    sequences = []
    seen_sequences = set()
    for t in ts:
        mRNA = db[g][t]['mRNA']
        if db[g][t].get("CDS"):
            for c in db[g][t]['CDS']:
                chrom, start, end, strand = c.chr, c.start, c.end, c.strand
                s = str(genome[chrom][start-1:end])
                s = Seq(s)
                if mRNA.strand == '-':
                    s = s.reverse_complement()
                s = str(s).upper()
                if s not in seen_sequences:
                    seen_sequences.add(s)
                    sequences.append(s)
        else:
            for e in db[g][t]['exon']:
                chrom, start, end, strand = e.chr, e.start, e.end, e.strand
                s = str(genome[chrom][start - 1:end])
                s = Seq(s)
                if mRNA.strand == '-':
                    s = s.reverse_complement()
                s = str(s).upper()
                if s not in seen_sequences:
                    seen_sequences.add(s)
                    sequences.append(s)
        gene_isoform[g] = sequences
    return gene_isoform

            
def get_sgrna(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # logging
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="sgRNA.log"
    )

    if args.pam:
        pam_table()
        return
    db, t2g = load_annotation_database(args.database)
    genome = genome_dict(args.genome)
    sgrna_seq_list = []
    # sgRNAs bind genomic DNA.  A source-aware SQLite catalogue prevents an
    # accidental transcriptome database from being used for off-target counts.
    kmer_bin_database = open_kmer_database(
        args.kmer_bin,
        k=args.length,
        expected_source='genome',
        jellyfish=getattr(args, 'jellyfish', None),
    )
    logging.info(f"Perform sgRNA design at the whole-genome scale")
    if not args.gene:
        for g in tqdm(db, ncols=80, total=len(db),
                      desc = 'sgRNA processing:', colour="GREEN"):
            gene = db[g]['gene']
            symbol = get_gene_symbol(db, g)
            if isinstance(gene, dict):
                print(gene, g)
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            # 
            ts = [i for i in db[g] if i not in {'gene', 'symbol'}]
            params = (g, ts, db, genome, args.style, args.output_format)
            gene_isoform = parse_cds_exon(params)
            test_transcripts = gene_isoform
            if len(test_transcripts) == 0:
                continue
            designer = CDSsgRNADesigner(g, test_transcripts,
                                        kmer_bin_database,
                                        pam_pattern=args.pam_pattern,
                                        sg_length=args.length)
            sgrna_list = designer.design(require_all=True)

            # fasta (defalut)
            if not sgrna_list:
                # no sgRNA
                continue
            for sgrna in sgrna_list:
                # sgrna is a dict
                sequence = sgrna['sequence']
                strand = sgrna['strand'] # not same to gene.strand
                gc = sgrna['gc_content']
                score = sgrna['score']
                targets = sgrna['targets']
                pam = sgrna['pam']
                kmer_freq = sgrna['kmer_freq']
                positions = sgrna['positions']

                if args.output_format == 'fasta':
                    desc='strand:%s transcript:%s PAM:%s start:%d end:%d length:%d score:%d kmer_freq:%d'%(gene.strand,
                                                                                                            "|".join(ts),
                                                                                                            pam,
                                                                                                            gene.start,
                                                                                                            gene.end,
                                                                                                            len(re.sub(r'[^ATCGatcg]', '', sequence)),
                                                                                                            score,
                                                                                                           kmer_freq)
                    sgRecord = SeqRecord(Seq(sequence), id=g, description=desc)
                    sgrna_seq_list.append(sgRecord)
                else:
                    it = [g, symbol, "|".join(ts), sequence, gc, strand, score, targets, pam, kmer_freq]
                    sgrna_seq_list.append(dict((_sgRNA_CSV_HEADER[i], it[i]) for i in range(len(_sgRNA_CSV_HEADER))))
        parse_output(args, sgrna_seq_list)
        logging.info(f"Finished sgRNA design at the whole-genome scale")
    else:
        logging.info(f"Perform sgRNA design ath the gene level")
        gene = db[args.gene]['gene']
        symbol = get_gene_symbol(db, args.gene)
        chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
        ts = [i for i in db[args.gene] if i not in {'gene', 'symbol'}]
        params = (args.gene, ts, db, genome, args.style, args.output_format)
        gene_isoform = parse_cds_exon(params)
        test_transcripts = gene_isoform
        # 设计靶向所有转录本的sgRNA
        designer = CDSsgRNADesigner(args.gene, test_transcripts,
                                    kmer_bin_database,
                                    pam_pattern=args.pam_pattern,
                                    sg_length=args.length)
        sgrna_list = designer.design(require_all=True)
        # print(args.gene, sgrna_list)
        # print("length:", len(sgrna_list))
        for sgrna in sgrna_list:
            sequence = sgrna['sequence']
            strand = sgrna['strand']  # not same to gene.strand
            gc = sgrna['gc_content']
            score = sgrna['score']
            targets = sgrna['targets']
            pam = sgrna['pam']
            kmer_freq = sgrna['kmer_freq']
            positions = sgrna['positions']
            if args.output_format == 'fasta':
                desc = 'strand:%s transcript:%s PAM:%s start:%d end:%d length:%d score:%d kmer_freq:%d' % (gene.strand,
                                                                                                              "|".join(ts),
                                                                                                              pam,
                                                                                                              gene.start,
                                                                                                              gene.end,
                                                                                                              len(re.sub(r'[^ATCGatcg]', '', sequence)),
                                                                                                              score,
                                                                                                           kmer_freq)
                sgRecord = SeqRecord(Seq(sequence), id=args.gene, description=desc)
                sgrna_seq_list.append(sgRecord)
            else:
                it = [args.gene, symbol, "|".join(ts), sequence, gc, strand, score, targets, pam, kmer_freq]
                sgrna_seq_list.append(dict((_sgRNA_CSV_HEADER[i], it[i]) for i in range(len(_sgRNA_CSV_HEADER))))
        parse_output(args, sgrna_seq_list)
        logging.info(f"Finished sgRNA design at the gene scale")

if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": "ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_002": "ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_003": "ATGAGATCCTGGCCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA"}
    designer = CDSsgRNADesigner("TP53", test_transcripts, kmer_bin=None)
    print(designer)
    # 寻找通用sgRNA
    universal = designer.design(require_all=True)
    print(f"找到 {len(universal)} 条通用sgRNA:")
    for sg in universal[:3]:
        print(f"Seq: {sg['sequence']} ({sg['strand']}链)")
        print(f"GC: {sg['gc_content']}, 评分: {sg['score']}")
        print(f"靶向转录本: {len(sg['positions'])}个\n")

    # 寻找任意sgRNA
    any_sg = designer.design(require_all=False)
    print(f"\n总候选sgRNA: {len(any_sg)} 条")
