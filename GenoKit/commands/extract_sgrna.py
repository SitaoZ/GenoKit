# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import parse_output, gff_feature_dict, gtf_feature_dict
from GenoKit.commands.design_sgrna import CDSsgRNADesigner

_CSV_HEADER = ['GeneID','TranscriptID','siRNA','GC','Strand']

def parse_transcript(params):
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    gene_isoform = defaultdict(lambda: defaultdict(dict))
    g, ts, db, genome, style, output_format = params
    transcript_seq_list = deque()
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
        seq = str(seq)
        gene_isoform[g][t] = seq
    return gene_isoform
        
            
def get_sgrna(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    #if args.style == 'GFF':
    #    db, t2g = gff_feature_dict(args.database, args.style)
    #else:
    #    db, t2g = gtf_feature_dict(args.database, args.style)
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    sirna_seq_list = []
    index = 0
    if not args.gene:
        for g in tqdm(db, ncols=80, total=len(db), desc = 'Gene processing:'):
            gene = db[g]['gene']
            if isinstance(gene, dict):
                print(gene, g)
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            # 
            ts = [i for i in db[g] if i != 'gene']
            params = (g, ts, db, genome, args.style, args.output_format)
            gene_isoform = parse_transcript(params)
            test_transcripts = gene_isoform[g]
            designer = SirnaDesigner(g, test_transcripts, length = 22)
            universal = designer.design(require_all=True)
            sirna_list = designer.get_top_designs(10)
            if args.output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                              gene.start,
                                                              gene.end,
                                                              len(seq))
                geneRecord = SeqRecord(seq, id=g, description=desc)
                gene_seq_list.append(geneRecord)
            else:
                for seq,gc in sirna_list:
                    it = [g, "|".join(ts), seq, gc, gene.strand]
                    sirna_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        parse_output(args, sirna_seq_list)
    else:
        gene = db[args.gene]['gene']
        chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
        ts =  [i for i in db[args.gene] if i != 'gene']
        params = (args.gene, ts, db, genome, args.style, args.output_format)
        gene_isoform = parse_transcript(params)
        test_transcripts = gene_isoform[args.gene]
        # 设计靶向所有转录本的siRNA
        designer = SirnaDesigner(args.gene, test_transcripts, length = 22)
        universal = designer.design(require_all=True)
        sirna_list = designer.get_top_designs(10)
        for seq,gc in sirna_list:
            it = [args.gene, "|".join(ts), seq, gc, gene.strand]
            sirna_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        parse_output(args, sirna_seq_list)
    

if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": "ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCAGGTTGTTCTTCTATGTGGCCCTTTGGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCTTTTGAAAGCAGTGATATAACTCTAGGTAAATGCTATGTCTACTAATTATAGTTTCTTAATTTTCATAGCTATATTATGAAAAGAGTAAATTGAAGAAATGGAAACTGCAAATTACACCAAGGTGACAGAATTTGTTCTCACTGGCCTATCCCAGACTCCAGAGGTCCAACTAGTCCTATTTGTTATATTTCTATCCTTCTATTTGTTCATCCTACCAGGAAATATCCTTATCATTTGCACCATCAGTCTAGACCCTCATCTGACCTCTCCTATGTATTTCCTGTTGGCTAATCTGGCCTTCCTTGATATTTGGTACTTAGCTAGGCTAGCTAGCTAGCTAGCTAGCGCGATCGATCGATCGATCGA",
        "NM_002": "ATCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        "NM_003": "ATGAGATCCTGGCCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGGCCCTTTGTCTCCTGTGGACAGGACACATGGATGCTGGAATCACCCAGAGCCCAAGACACAAGGTCACAGAGACAGGAACACCAGTGACTCTGAGATGTCACCAGACTGAGAACCACCGCTATATGTACTGGTATCGACAAGACCCGGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGTTAAAGATACTGACAAAGGAGAAGTCTCAGATGGCTATAGTGTCTCTAGATCAAAGACAGAGGATTTCCTCCTCACTCTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTC"}
    designer = CDSsgRNADesigner("TP53", test_transcripts)

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
