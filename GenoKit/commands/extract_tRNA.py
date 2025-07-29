# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.database.database import create_db



def get_tRNA_gb(args):
    tRNA = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'tRNA': # CDS promoter UTR 
                tRNA_seq = ''
                for part in feature.location.parts:
                    tRNA_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                product = feature.qualifiers['product'][0] if 'product' in feature.qualifiers else 'Null'
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                tRNA_seq_record = SeqRecord(tRNA_seq, 
                                            id='gene:%s'%(gene_id),
                                            description='%s strand %s length %d'%(product,feature.strand, len(tRNA_seq))
                                           )
                tRNA.append(tRNA_seq_record)
                #break 
    if args.print and args.format == 'dna':
        SeqIO.write(tRNA, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(tRNA, args.output, 'fasta')
