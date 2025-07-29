# -*- coding: utf-8 -*-
import sys
import time
import gffutils
import pandas as pd 
import _pickle as cPickle
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.utils.util import mRNA_type, parse_output
from GenoKit.database.database import create_db, genome_dict

_CSV_HEADER = ['TranscriptID','Chrom','Start','End','Strand','Exon', 'Exon_length', 'Exon_index']

def sub_exon(params):
    """
    parameters:
        params:
    return:
        exon_seq_list, a deque class 
    """
    g, t, tdb, genome, style, output_format = params
    exon_seq_list = deque()
    exon_index = 0
    
    if output_format == 'gff':
        # gff
        for line in tdb['exon']:
            exon_seq_list.append("\t".join(map(str, list(line))))
    else:
        mRNA = tdb['mRNA']
        for e in tdb['exon']:
            seq = genome[e.chr][e.start-1 : e.end].seq
            seq = Seq(seq)
            if e.strand == '-':
                seq = seq.reverse_complement()
            exon_index += 1
            if output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s exon %s start:%d end:%d length=%d'%(e.strand,
                                                                  exon_index,
                                                                  e.start,
                                                                  e.end,
                                                                  len(seq))
                exonRecord = SeqRecord(seq, id=t, description=desc)
                exon_seq_list.append(exonRecord)
            else:
                # csv
                it = [t, e.chr, e.start, e.end, e.strand, seq, len(seq), exon_index]
                exon_seq_list.append(dict(zip(_CSV_HEADER, it)))
    return exon_seq_list    



def get_exon(args):
    '''
    parameters:
     
    '''
    starttime = time.time()
    #if args.style == 'GFF':
    #    db, t2g = gff_feature_dict(args.database, args.style)
    #else:
    #    db, t2g = gtf_feature_dict(args.database, args.style)
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    endtime = time.time()
    print(endtime - starttime)
    starttime = time.time()
    genome = genome_dict(args.genome)
    endtime = time.time()
    print(endtime - starttime)
    
    if not args.transcript:
        param_list = [(g, t, db[g][t], genome, args.style, args.output_format) for g in db for t in db[g] if t != 'gene']
        exon_seq_list = deque()
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='Exon processing:'):
            exon_seq_list.append(sub_exon(para))
        
        exon_seq_list = [d for de in exon_seq_list if de != None for d in de]
        starttime = time.time()
        #if args.output_format == 'fasta':
        #    with open(args.output, 'w') as handle:
        #        for record in exon_seq_list:
        #            handle.write(">" + record.id + " " + record.description + "\n" + str(record.seq) + "\n")
        #else:
        #    parse_output(args, exon_seq_list)
        parse_output(args, exon_seq_list)
        endtime = time.time()
        print(endtime - starttime)
    else:
        # g, t, tdb, genome, style, output_format = params
        g = t2g[args.transcript]
        t = args.transcript
        para = (g, t, db[g][t], args.genome, args.style, args.output_format)
        exon_seq_list = sub_exon(para)
        exon_seq_list = [d for de in exon_seq_list if de != None for d in de]
        # output
        parse_output(args, exon_seq_list)

def get_exon_gb(args):
    exons = []
    for record in create_db(args.genbank):
        index = 1
        print(type(record))
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                exon_seq = feature.extract(record.seq)
                # feature.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                exon_seq_record = SeqRecord(exon_seq,id='%s exon %d'%(record.id, index),
                                 description='strand %s length %d'%(feature.strand, len(exon_seq)))
                exons.append(exon_seq_record)
                index += 1 
    if args.print:
        SeqIO.write(exons, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(exons, args.output, 'fasta')
 
