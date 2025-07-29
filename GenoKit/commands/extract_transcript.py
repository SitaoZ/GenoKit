# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
import _pickle as cPickle
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from GenoKit.database.database import create_db, genome_dict
from GenoKit.utils.util import add_stop_codon, mRNA_type, seq_upper_lower, parse_output

'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''
def stop_codon(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s


def anchor_CDS(g, t, db, genome, style):
    '''
    parameters: 
     db : database create by gffutils
     transcript: gffutils transcript feature type
     style: database file type
    return:
     start_codon position 
     stop_codon positiopn 
     cds length 
    '''
    mRNA = db[g][t]['mRNA']
    cds = ''
    cds_len = 0
    if db[g][t].get('CDS'):
        for c in db[g][t].get('CDS'):
            chrom, start, end, strand = c.chr, c.start, c.end, c.strand
            cds += str(genome[chrom][start-1:end])
        cds = Seq(cds)
        if mRNA.strand == '-':
            cds = cds.reverse_complement()
        cds_len = len(cds) # add stop codon length
    utr5 = ''
    utr5_len = 0
    if db[g][t].get('five_prime_UTR'):
        for u in db[g][t].get('five_prime_UTR'):
            chrom, start, end, strand = u.chr, u.start, u.end, u.strand
            utr5 += str(genome[chrom][start-1:end])
        utr5_len = len(utr5)
    utr3 = ''
    utr3_len = 0
    if db[g][t].get('three_prime_UTR'):
        for u in db[g][t].get('five_prime_UTR'):
            chrom, start, end, strand = u.chr, u.start, u.end, u.strand
            utr3 += str(genome[chrom][start-1:end])
        utr3_len = len(utr3)
          
    return utr5_len, utr3_len, cds_len


_CSV_HEADER = ['TranscriptID', 'Chrom', 'Start_genome', \
               'End_genome','Start_transcript', \
               'End_transcript','Strand','CDS']

def sub_transcript(params):
    """
    parameters:
        params: argparse
    return:
        transcript_seq_list, a deque
    """
    g, t, db, genome, style, upper, output_format = params
    transcript_seq_list = deque()
    utr5_len, utr3_len, cds_len = anchor_CDS(g, t, db, genome, style)
    atg2firstexon = 0
    seq = ''
    mRNA = db[g][t]['mRNA']
    for e in db[g][t]['exon']:
        chrom, start, end, strand = e.chr, e.start, e.end, e.strand
        s = str(genome[chrom][start-1:end])
        seq += s
    seq = Seq(seq)
    if mRNA.strand == '-':
        seq = seq.reverse_complement()
    if upper:
        seq = seq_upper_lower(seq, utr5_len , utr5_len + cds_len)
    # position in transcript 1 based 
    if db[g][t].get('CDS'):
        cds_start_transcript = utr5_len + 1
    else:
        cds_start_transcript = 0
    cds_end_transcript = utr5_len + cds_len
    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(mRNA.strand, mRNA.start, mRNA.end, len(seq),
                                                              cds_start_transcript,
                                                              cds_end_transcript)
        transcriptRecord = SeqRecord(seq, id=t, description=desc)
        transcript_seq_list.append(transcriptRecord)
    else:
        it = [t, mRNA.chr, mRNA.start, mRNA.end,\
              cds_start_transcript, cds_end_transcript, mRNA.strand, seq]
        transcript_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return transcript_seq_list

    
def get_transcript(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    genome = genome_dict(args.genome)
    # feature_types = db.featuretypes()
    # assert GTF or GFF
    #if args.rna_feature == 'mrna':
    #    mRNA_str = mRNA_type(args.style)
    #else:
    #    mRNA_str = tuple(x for x in list(feature_types) if 'RNA' in x or 'transcript' in x)
    if not args.transcript:
        param_list = [(g, t, db, genome, args.style, args.upper, args.output_format) for g in db for t in db[g] if t != 'gene']
        transcript_seq_list = deque()
        for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='Transcript processing:'):
            transcript_seq_list.append(sub_transcript(para))
        transcript_seq_list = [d for de in transcript_seq_list if de != None for d in de]
        parse_output(args, transcript_seq_list)
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                param = (t, args.database, args.genome, args.style, args.upper, args.output_format)
                transcript_seq_list = sub_transcript(param)
                parse_output(args, transcript_seq_list)
                break 


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_transcript_gb(args):
    '''
    '''
    transcript = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                transcript_seq = ''
                #for part in feature.location.parts:
                #    cds_seq += part.extract(record.seq)
                transcript_seq = feature.location.extract(record).seq
                # part.strand 会将FeatureLocation -1的反向互补
                # geneid or locus_tag 
                if 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                else:
                    gene_id = "Null"
                # protein id  
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                else:
                    protein_id = 'Null'

                transcript_seq_record = SeqRecord(cds_seq,id='gene:%s protein:%s'%(gene_id, protein_id),
                                 description='strand %s length %d'%(feature.strand, len(cds_seq)))
                transcript.append(cds_seq_record)
                break 
        break     
    if args.print and args.format == 'dna':
        SeqIO.write(cds, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(proteins, sys.stdout, 'fasta')
    elif args.output and args.format == 'dna':
        SeqIO.write(cds, args.output, 'fasta')
    elif args.output and args.format == 'protein':
        SeqIO.write(proteins, args.output, 'fasta')
