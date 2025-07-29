# -*- coding: utf-8 -*-
import io
import os
import sys
import time
import numba
import gffutils
import threading
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
import _pickle as cPickle
from multiprocessing import Pool, Manager
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque, Counter
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import add_stop_codon, mRNA_type, parse_output, gff_feature_dict, gtf_feature_dict



def get_stat(args):
    '''
    parameters:
        args: parse from argparse
    return:
        a genome database statistics
    '''
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    chrom_size = len(genome.keys())
    genome_size = len(genome)
    gene_nums = len(db)
    transcript_nums = len([t for g in db for t in db[g] if t != 'gene'])
    coding_gene_nums = len(set([g for g in db for t in db[g] if t != 'gene' and db[g][t].get('CDS')]))
    coding_transcript_nums = len(set([t for g in db for t in db[g] if t != 'gene' and db[g][t].get('CDS')]))
    utr5_transcript_nums = len(set([t for g in db for t in db[g] if t != 'gene' and db[g][t].get('five_prime_UTR')])) 
    utr3_transcript_nums = len(set([t for g in db for t in db[g] if t != 'gene' and db[g][t].get('three_prime_UTR')]))
    multiple_exon_transcript_nums = len(set([t for g in db for t in db[g] if t != 'gene' and db[g][t].get('exon') and len(db[g][t]['exon']) >=2]))
    print('Genome      size:', chrom_size)
    print('Chromosome  size:', genome_size)
    print('Gene       count:', gene_nums)
    print('Transcript count:', transcript_nums)
    print('Gene with CDS:', coding_gene_nums)
    print('Transcript with CDS:', coding_transcript_nums)
    print('Transcript with UTR5:', utr5_transcript_nums)
    print('Transcript with UTR3:', utr3_transcript_nums)
    print('Transcript with multiple exon:', multiple_exon_transcript_nums)
