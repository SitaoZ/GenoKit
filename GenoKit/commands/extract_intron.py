# -*- coding: utf-8 -*-
import sys
import time
import gffutils
import itertools
import pandas as pd 
import _pickle as cPickle
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import genome_dict
from GenoKit.utils.util import utr3_type, utr5_type, mRNA_type, parse_output

_CSV_HEADER = ['TranscriptID', 'Intron', 'Chrom','Start_genome', \
               'End_genome', 'Strand', 'Sequence']

def sub_intron(params):
    '''
    parameters:
        params: argparse
    return:
        intron_seq_list, a deque
    '''
    # g, t, db[g][t], genome, args.style, args.output_format
    g, t, tdb, genome, style, output_format = params
    mRNA = tdb['mRNA']
    intron_seq_list = deque()
    exons = list(itertools.chain(*[[e.start,e.end] for e in tdb['exon']])) # exon position include start and end
    if len(exons)/2 == 1:
        # only one exon
        return None
    else:
        left = exons[1:-1]
        intron_index = 0
        for i in range(0, len(left), 2):
            s = left[i]   # intron start
            e = left[i+1] # intron end, shoukd be minus 1 
            if s+1 != e:
                # exon 与exon 不相连，就是intron的位置
                # 注意索引的点，
                intron_index += 1
                intron_seq= genome[mRNA.chr][s:e-1]
                if mRNA.strand == '-':
                    intron_seq = intron_seq.reverse
                    intron_seq = intron_seq.complement
                    intron_rank = int(len(exons)/2 - intron_index)
                else:
                    intron_rank = intron_index
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    intron_gff = [mRNA.chr, 'GenoKit', 'intron', s+1, 
                                  e-1, '.', mRNA.strand, '.', f'ID={t};Name=intron{intron_rank};']
                    intron_gff_line = list(map(str, intron_gff))
                    intron_seq_list.append("\t".join(intron_gff_line))
                    continue
                elif output_format == 'fasta':
                    intron_seq= Seq(str(intron_seq))
                    # 1 based position
                    intronRecord = SeqRecord(intron_seq, id=t, 
                                             description='strand %s intron%d start %d end %d length=%d exons=%d'%(
                                             mRNA.strand, intron_rank, s+1, e-1, len(intron_seq), len(exons)/2)
                                            )
                    intron_seq_list.append(intronRecord)
                else:
                    intron_seq= str(intron_seq)
                    it = [t, f'intron{intron_rank}', mRNA.chr, s+1, e-1, mRNA.strand, intron_seq]
                    intron_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        return intron_seq_list

def get_intron(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
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
    
    intron_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Intron']) # header 
    if not args.transcript:
        param_list = [(g, t, db[g][t], genome, args.style, args.output_format) for g in db for t in db[g] if t != 'gene']
        intron_seq_list = deque()
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='Intron processing:'):
            intron_seq_list.append(sub_intron(para))

        intron_seq_list = [d for de in intron_seq_list if de != None for d in de]
        starttime = time.time()
        parse_output(args, intron_seq_list)
        endtime = time.time()
        print(endtime - starttime)
    else:
        g = t2g[args.transcript]
        t = args.transcript
        para = (g, t, db[g][t], genome, args.style, args.output_format) 
        intron_seq_list = sub_intron(para)
        intron_seq_list = [d for de in intron_seq_list if de != None for d in de]
        parse_output(args, intron_seq_list)
