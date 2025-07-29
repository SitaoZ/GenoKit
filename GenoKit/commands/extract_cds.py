# -*- coding: utf-8 -*-
import sys
import time
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
from multiprocessing import Pool, Process, Queue
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import genome_dict
from GenoKit.utils.util import add_stop_codon, mRNA_type, parse_output, gff_feature_dict, gtf_feature_dict

#import heartrate
#heartrate.trace(browser=True)
'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''

_CSV_HEADER = ['TranscriptID','Chrom','Start','End','Strand','CDS']

def sub_cds(params):
    g, t, db, genome, style, output_format = params
    cds_seq_list = deque()
    seq = ''
    cds_start_transcript = 0
    cds_end_transcript = 0
    cds_gff_lines = deque()
    mRNA = db[g][t]['mRNA']
    if db[g][t].get('CDS'):
        for c in db[g][t]['CDS']:
            cds_gff_lines.append(c)
            chrom, start, end, strand = c.chr, c.start, c.end, c.strand
            s = str(genome[chrom][start-1:end])
            seq += s
            if not cds_start_transcript:
                cds_start_transcript = abs(start - mRNA.start)
            cds_end_transcript += len(s)
        cds_end_transcript = cds_end_transcript + cds_start_transcript
    else:
        return
    seq = Seq(seq)
    if mRNA.strand == '-':
        seq = seq.reverse_complement()
    if output_format == 'gff':
        # gff
        for line in cds_gff_lines:
            cds_seq_list.append(line)
    elif output_format == 'fasta':
        # fasta (defalut)
        desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand,
                                                              mRNA.start,
                                                              mRNA.end,
                                                              len(seq),
                                                              cds_start_transcript,
                                                              cds_end_transcript)
        cdsRecord = SeqRecord(seq, id=t.replace('transcript:',''), description=desc)
        cds_seq_list.append(cdsRecord)
    else:
        # csv
        it = [t, mRNA.chr, mRNA.start, mRNA.end, mRNA.strand, seq]
        cds_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return cds_seq_list


def sub_cds2(params):
    g, t, tdb, genome, style, output_format = params
    # genome = genome_dict(genome)
    cds_seq_list = deque()
    if output_format == 'gff':
        # gff
        cds_gff_lines = (c for c in tdb['CDS'])
        for line in cds_gff_lines:
            cds_seq_list.append("\t".join(map(str, list(line))))
    else:
        mRNA = tdb['mRNA']
        cds_start_position = tdb['CDS'][0].start - mRNA.start if mRNA.start == '+' else mRNA.end - tdb['CDS'][-1].end
        seq = [genome[c.chr][c.start-1 : c.end].seq for c in tdb['CDS']]
        seq = ''.join(seq)
        cds_end_position = len(seq) + cds_start_position
        seq = Seq(seq)
        if mRNA.strand == '-':
            seq = seq.reverse_complement()
        
        if output_format == 'fasta':
            # fasta (defalut)
            desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand,
                                                                  mRNA.start,
                                                                  mRNA.end,
                                                                  len(seq),
                                                                  cds_start_position,
                                                                  cds_end_position)
            cdsRecord = SeqRecord(seq, id=t.replace('transcript:',''), description=desc)
            cds_seq_list.append(cdsRecord)
        else:
            # csv it = [t, mRNA.chr, mRNA.start, mRNA.end, mRNA.strand, seq]
            it = {'TranscriptID': t,
                  'Chrom': mRNA.chr,
                  'Start': mRNA.start,
                  'End' : mRNA.end,
                  'Strand' : mRNA.strand,
                  'CDS' : seq}
            cds_seq_list.append(it)
    return cds_seq_list

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def split_list(l: list, parts: int) -> list:
    """Takes a list as input, and splits it into "parts" number of sub-lists,
    which are then inserted as elements in the returned meta-list.
    
    The function will try to make the sub-lists as equal in length as 
    possible, so running
    split_list( [1, 2, 3, 4, 5, 6], 4 ) 
    will return 
    [ [1, 2], [3, 4], [5], [6] ]
    
    I will also make sure that the list isn't split into more parts than
    there are elements in the list, so
    split_list( [a, b], 6 ) 
    will return 
    [ [a], [b] ]
    """
    n = min(parts, max(len(l),1))
    k, m = divmod(len(l), n)
    return [l[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]




def multi_mediator(param_list):
    cds_seq_list = deque()
    # str to database
    genome = genome_dict(param_list[0][3])
    for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='CDS processing', colour="red"):
        para = list(para)
        para[3] = genome # add genome pyfaidx 
        cds_seq_list.append(sub_cds2(para))
    return cds_seq_list


def multi_mediator2(param_list, queue):
    cds_seq_list = deque()
    # str to database
    genome = genome_dict(param_list[0][3])
    for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='CDS processing', colour="red"):
        para = list(para)
        para[3] = genome # add genome pyfaidx 
        cds_seq_list.append(sub_cds2(para))
    queue.put(cds_seq_list)
    return cds_seq_list

def get_cds(args):
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
    
    if not args.transcript:
        # single cpu 
        starttime = time.time()
        genome = genome_dict(args.genome)
        endtime = time.time()
        print(endtime - starttime)
        param_list = [(g, t, db[g][t], genome, args.style, args.output_format) for g in db for t in db[g] if t != 'gene' and db[g][t].get('CDS')]
        cds_seq_list = deque()
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='CDS processing', colour="GREEN"):
            cds_seq_list.append(sub_cds2(para))
        
        cds_seq_list = [d for de in cds_seq_list if de != None for d in de]
        ## multiprocessing
        """
        param_list = [(g, t, db[g][t], args.genome, args.style, args.output_format) for g in db for t in db[g] if t != 'gene' and db[g][t].get('CDS')]
        list_parts = list(np.array_split(param_list, args.process))
        with Pool(processes=args.process) as p:
            cds_seq_list = p.map(multi_mediator, list_parts)
        cds_seq_list = [d[0] for de in cds_seq_list if de != None for d in de]
        """

        """
        # method 2
        queue = Queue()
        for i in range(args.process):
            p = Process(target=multi_mediator2, args=(list_parts[i], queue))
            p.start()
            p.join()
        cds_seq_list = queue.get() 
        """
        starttime = time.time()
        if args.output_format == 'fasta':
            with open(args.output, 'w') as handle:
                for record in cds_seq_list:
                    handle.write(">" + record.id + " " + record.description + "\n" + str(record.seq) + "\n")
        else:
            parse_output(args, cds_seq_list)
        endtime = time.time()
        #print(endtime - starttime)
    else:
        # only one transcript
        if not t2g.get(args.transcript):
            sys.exit(1)
            print('transcript id not exist in database, please check.')
        param = (t2g[args.transcript], args.transcript, db, genome, args.style, args.output_format)
        cds_seq_list = sub_cds(param)
        parse_output(args, cds_seq_list)
