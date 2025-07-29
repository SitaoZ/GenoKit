# -*- coding: utf-8 -*-
import re
import sys
import pickle
import gffutils
import itertools
import pandas as pd 
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import parse_output, gff_feature_dict
from GenoKit.utils.util import record

_CSV_HEADER = ['IGR_ID', 'Chrom', 'Strand', 'Start_genome', \
               'End_genome', 'Sequence']

def sub_IGR(params):
    """
    parameters:
        params:
    return:
        igr_seq_list, a deque class 
    """
    chrom, database, genome, style, output_format, genes = params
    if style == 'GFF':
        g_pat = re.compile(r'ID=(\S+?);')
    else:
        g_pat = re.compile(r'gene_id "(\S+?)";')
    igr_seq_list = deque()
    gene_positions = list(itertools.chain(*[[g.start,g.end] for g in genes]))
    # gene_ids = [g.attribute for g in genes]
    gene_ids = []
    for g in genes:
        match = g_pat.search(g.attribute)
        gene_ids.append(match.group(1)) 
    strands = [g.strand for g in genes]
    if len(gene_positions)/2 == 1:
        # only one gene
        return None
    else:
        left = gene_positions[1:-1] # index 1
        for i in range(0, len(left), 2):
            s = left[i]   # igr start
            e = left[i+1] # igr end, shoukd be minus 1 
            igr_id = gene_ids[int(i/2)]+"_"+gene_ids[int(i/2)+1]
            igr_strand = strands[int(i/2)]+"|"+strands[int(i/2)+1]
            if s+1 != e:
                # gene 与gene 不相连，就是igr的位置
                # 注意索引的点，
                if s > e-1:
                    # gene overlap
                    continue 
                igr_seq= genome[chrom][s:e-1]
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    igr_gff = [chrom, 'GenoKit', 'igr', s+1, e-1, '.', igr_strand, '.', f'ID={igr_id};']
                    igr_gff_line = list(map(str, igr_gff))
                    igr_seq_list.append("\t".join(igr_gff_line))
                    continue
                elif output_format == 'fasta':
                    igr_seq= Seq(str(igr_seq))
                    # 1 based position
                    igrRecord = SeqRecord(igr_seq, id=f'{igr_id}',
                                             description='strand %s start %d end %d length=%d'%(
                                             igr_strand, s+1, e-1, len(igr_seq))
                                            )
                    igr_seq_list.append(igrRecord)
                else:
                    igr_seq= str(igr_seq)
                    it = [igr_id, chrom, igr_strand, s+1, e-1, igr_seq]
                    igr_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        return igr_seq_list

def get_igr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # GTF noncoding are annotated transcript,while GFF noncoding ared annotated gene
    #if args.style == 'GFF':
    #    db, t2g = gff_feature_dict(args.database, args.style)
    #else:
    #    db, t2g = gtf_feature_dict(args.database, args.style)
    with open(args.database, 'rb') as f:
        db = pickle.load(f)
    genome = genome_dict(args.genome)
    chrom_split = dict()
    for g in db:
        gene = db[g]['gene']
        chrom = gene.chr
        if chrom_split.get(chrom):
            chrom_split[chrom].append(gene)
        else:
            chrom_split[chrom] = [gene]
    param_list = [(chrom, args.database, genome, args.style, args.output_format, chrom_split[chrom]) for chrom in chrom_split.keys()]
    igr_seq_list = deque()
    for param in tqdm(param_list, total=len(param_list), ncols=80, desc='IGR processing:'):
        igr_seq_list.append(sub_IGR(param))
    igr_seq_list = [d for de in igr_seq_list if de != None for d in de]
    parse_output(args, igr_seq_list)

def get_IGR_gb(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    gene = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'gene': # CDS promoter UTR 
                gene_seq = ''
                for part in feature.location.parts:
                    gene_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                gene_seq_record = SeqRecord(gene_seq, id='gene:%s'%(gene_id), 
                                  description='strand %s length %d'%(feature.strand, len(gene_seq)))
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
