# -*- coding: utf-8 -*-
import os 
from Bio import SeqIO
from pyfaidx import Fasta

def genome_dict(genome_fasta_path):
    '''
    parameters:
     genome_fasta_path: genome reference of organism
    return:
     genome fasta dict
    '''
    # genome = Fasta(genome_fasta_path)
    genome = Fasta(genome_fasta_path, rebuild=False, build_index=False)
    return genome

def create_db(genbank_path):
    if not os.path.exists(genbank_path):
        print('args.genbank path not exist; please check.')
        sys.exit(1)
    for record in SeqIO.parse(genbank_path,'gb'):
        yield record
