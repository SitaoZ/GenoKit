# -*- coding: utf-8 -*-
'''
filename: database.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
import sys
from Bio import SeqIO
from pyfaidx import Fasta, Faidx

def genome_dict(genome_fasta_path):
    '''
    parameters:
     genome_fasta_path: genome reference of organism
    return:
     genome fasta dict
    '''

    # Fasta faidx
    if not os.path.exists(genome_fasta_path + ".fai"):
        # Faidx automatically creates a missing .fai index.  Do not invoke
        # build_index() a second time, which would reread the full FASTA.
        Faidx(genome_fasta_path)
    # genome = Fasta(genome_fasta_path)
    # genome file format only support fasta an BGZIP file .
    genome = Fasta(genome_fasta_path, rebuild=False, build_index=False)
    return genome

def gb_iter(genbank_path):
    if not os.path.exists(genbank_path):
        print('args.genbank path not exist; please check.')
        sys.exit(1)
    for record in SeqIO.parse(genbank_path,'gb'):
        yield record
