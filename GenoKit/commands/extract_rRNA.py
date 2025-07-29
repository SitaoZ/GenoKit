# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.database.database import create_db



def get_rRNA_gb(args):
    rRNA = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'rRNA': # CDS promoter UTR 
                rRNA_seq = ''
                for part in feature.location.parts:
                    rRNA_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                product = feature.qualifiers['product'][0] if 'product' in feature.qualifiers else 'Null'
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                rRNA_seq_record = SeqRecord(rRNA_seq, 
                                            id='gene:%s'%(gene_id), 
                                            description='%s strand %s length %d'%(product,feature.strand, len(rRNA_seq))
                                           )
                rRNA.append(rRNA_seq_record)
                #break 
    if args.print and args.format == 'dna':
        SeqIO.write(rRNA, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(rRNA, args.output, 'fasta')
