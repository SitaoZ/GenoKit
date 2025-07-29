# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import parse_output, gff_feature_dict, gtf_feature_dict

_CSV_HEADER = ['GeneID','Chrom','Start','End','Strand','Seq']

def get_gene(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    # gene_seq = pd.DataFrame(columns=_CSV_HEADER)
    if args.style == 'GFF':
        db, t2g = gff_feature_dict(args.database, args.style)
    else:
        db, t2g = gtf_feature_dict(args.database, args.style)
    genome = genome_dict(args.genome)
    gene_seq_list = []
    index = 0
    if not args.gene:
        for g in tqdm(db, ncols=80, total=len(db), desc = 'Gene processing:'):
            gene = db[g]['gene']
            if isinstance(gene, dict):
                print(gene, g)
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            seq = str(genome[chrom][start-1:end])
            seq = Seq(seq)
            if gene.strand == '-':
                seq = seq.reverse_complement()
            it = [g, gene.chr, gene.start, gene.end, gene.strand, seq]
            if args.output_format == 'gff':
                gene_seq_list.append('\t'.join(it))
            elif args.output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                              gene.start,
                                                              gene.end,
                                                              len(seq))
                geneRecord = SeqRecord(seq, id=g, description=desc)
                gene_seq_list.append(geneRecord)
            else:
                gene_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        parse_output(args, gene_seq_list)
    else:
        gene = db[args.gene]['gene']
        chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
        seq = str(genome[chrom][start-1:end])
        seq = Seq(seq)
        if gene.strand == '-':
            seq = seq.reverse_complement()
        it = [args.gene, gene.chr, gene.start, gene.end, gene.strand, seq]
        if args.output_format == 'gff':
            gene_seq_list.append('\t'.join(it))
        elif args.output_format == 'fasta':
            # fasta (defalut)
            desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                          gene.start,
                                                          gene.end,
                                                          len(seq))
            geneRecord = SeqRecord(seq, id=g, description=desc)
            gene_seq_list.append(geneRecord)
        else:
            gene_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        parse_output(args, gene_seq_list)

def get_gene_gb(args):
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
                                            description='strand %s length %d'%(
                                            feature.strand, len(gene_seq))
                                           )
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
