#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: fasta_split.py
date: 2025/12/3 12:05
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
from tqdm import tqdm
from Bio import SeqIO
from pyfaidx import Fasta
from GenoKit.utils.util import open_auto


def get_fasta_bio(args):
    '''

    :param args:
    args.fasta: fasta file path
    args.input: input fasta ID
    args.output_dir: output dir
    args.prefix : output file prefix
    :return:
    '''
    seq_dict = dict()
    id_set = set()

    for record in tqdm(SeqIO.parse(args.fasta, 'fasta'), ncols=80, total = 100,
                       desc='Fasta processing:', colour="GREEN"):
        seq_dict[record.id] = str(record.seq)
        id_set.add(record.id)
    output = os.path.join(args.output_dir, f"{args.prefix}.{args.style}")
    with open(args.input, 'r') as f, open(output, 'w') as out:
        for line in f:
            line = line.strip()
            if line in id_set:
                out.write(">"+line+'\n')
                out.write(seq_dict[line]+'\n')

def fasta_count(filepath):
    fasta = Fasta(filepath)
    count = len(fasta.keys())
    return count

def count_lines(filepath):
    """快速获取文件行数"""
    with open_auto(filepath, mode='rb') as f:
        return sum(1 for _ in f)


def get_fasta(args):
    '''
    fasta 文件运行很快，fasta.gz运行很慢，主要在于f"{record}\n"
    :param args:
    args.fasta: fasta file path
    args.input: input fasta ID
    args.output_dir: output dir
    args.prefix : output file prefix
    :return:
    '''
    fasta = Fasta(args.fasta)
    id_set = set(fasta.keys())
    line_count = count_lines(args.input)
    with open_auto(args.input, 'r') as f, open(args.output, 'w') as out:
        for line in tqdm(f, ncols=80, total = line_count,
                       desc='Fasta processing:', colour="GREEN"):
            fasta_id = line.strip()
            if fasta_id in id_set:
                record = fasta[fasta_id]
                # sequence = record[:]
                out.write(">"+fasta_id+" "+ str(len(record))+'\n')
                out.write(f"{record}\n")