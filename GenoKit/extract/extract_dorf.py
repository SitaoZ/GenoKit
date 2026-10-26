# -*- coding: utf-8 -*-
'''
filename: extract_dorf.py
date: 2025/11/29 15:28
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
from tqdm import tqdm 
from Bio.Seq import Seq
from typing import Tuple
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from GenoKit.extract.extract_uorf import GFF
from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_transcript_models_with_symbols,
    load_annotation_database, transcript_count,
)
from GenoKit.utils.visualize import visual_transcript, visual_transcript_genome
from GenoKit.utils.util import parse_output, stop_codon_seq, add_stop_codon, stream_parse_output
from GenoKit.utils.gtf_gff_tool import gff_attr_to_gtf

# table header 
_dORF_csv_header = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Chrom', 'Strand',
    'CDS Interval', 'dORF Start', 'dORF End',
    'dORF Type', 'dORF Length', 'dORF'
]



class dORF(object):
    """dORF finder using object-oriented
    Attributes:
        transcript_id (str): transcript id
        length (int): orf length threshold
        chrom (str): chromosome id
        strand (str): transcript strand 
        matural_transcript (str): matural transcript sequence
        coding_sequence (str): conding sequence of specifity transcript
    """
    
    def __init__(self,
                 gene_id:str,
                 transcript_id: str,
                 length: int,
                 chrom: str,
                 strand: str,
                 matural_transcript:str,
                 coding_sequence:str):
        self.gene_id= gene_id
        self.transcript_id = transcript_id
        self.length = length
        self.chrom = chrom
        self.strand = strand 
        self.mt = matural_transcript
        self.cds = coding_sequence
        self.START_CODON = 'ATG'
        self.STOP_CODON_LIST = set(['TAA', 'TAG', 'TGA'])

    def transcript_location(self):
        """Transcript start and end
        1D list regradless of strand  
        """
        return [0, len(self.mt)]

    def cds_location_transcript(self):
        """Extract the cds start and end in transcript
        1D list regardless of strand 
        """
        start_codon_position = self.mt.index(self.cds)
        cds_start = start_codon_position + 1
        cds_end = start_codon_position + len(self.cds)
        return [cds_start, cds_end]

    def dorf_location_transcript(self):
        """Extract dorfs start and end in transcript
        2D list regardless strand 
        """
        d2 = []
        dORF_dict = self.dorf_parse()
        for dorf_type in dORF_dict.keys():
            for it in dORF_dict[dorf_type]:
                if len(it[9]) < self.length:
                    continue
                dorf_start, dorf_end = it[5], it[6]
                d2.append([dorf_start, dorf_end])
        return d2

    def find_start_codons(self, utr3, utr3_len):
        """find the start codon positions"""
        positions = []
        start_idx = 0
        while start_idx <= utr3_len - 3:
            pos = utr3.find(self.START_CODON, start_idx, utr3_len)
            if pos == -1:
                break
            positions.append(pos)
            start_idx = pos + 1  # 移动到下一个位置继续查找
        return positions

    def dorf_parse(self):
        """ dORF parse"""
        dORF_dict = defaultdict(list)
        mt = self.mt
        cds = self.cds
        # start_codon means the first base position in matural_transcript
        if not cds:
            return dORF_dict
        # print(self.transcript_id,"mt",mt)
        # print(self.transcript_id, "cds",cds)
        start_codon = mt.index(self.cds)
        # stop_codon means the last base position in matural transcript
        cds_len = len(cds)
        stop_codon = start_codon + cds_len
        cds_intervel = str(start_codon) + '-' + str(stop_codon)
        mt_len = len(mt)
        # utr5 = mt[:start_codon]
        # utr5_len = len(utr5)
        utr3 = mt[stop_codon:]
        utr3_len = len(utr3)
        find_start_codons_pos = self.find_start_codons(utr3, utr3_len)
        if len(find_start_codons_pos) == 0:
            return dORF_dict
        # method1
        # for i in range(0, utr3_len-3):
        # method2
        for i in find_start_codons_pos:
            # start codon find 
            # if utr3[i:i+3] == dORF.START_CODON:
            for j in range(i+3, utr3_len, 3):
                codon = utr3[j:j+3]
                # stop codon find
                if codon in self.STOP_CODON_LIST:
                    dorf = utr3[i:j+3] # str for polars write
                    out1 = [self.gene_id, self.transcript_id, self.chrom, self.strand,
                            cds_intervel, stop_codon+i+1, stop_codon+j+3,
                            'dORF', len(dorf), dorf]
                    if not dORF_dict.get('dORF'):
                        dORF_dict['dORF'] = [out1]
                    else:
                        dORF_dict['dORF'].append(out1)
                    break
        return dORF_dict
     
def dorf_procedure_oriented(transcript_id, chrom, strand, matural_transcript, coding_sequence):
    """
    Parameters:
        transcript_id (str):      transcript id
        chrom (str): chromsome id
        strand (str): strand
        matural_transcript: a Seq type (Biopython) from mature transcript without intron
        coding_sequence:    a Seq type (Biopython) from coding sequence start with ATG , 
                         end with TAA, TGA, TAG
    Return:
        down stream open reading frame in utr3
    """
    dORF_dict = defaultdict(list)
    stop_codon_list = ['TAA','TAG','TGA']
    # start_codon means the first base position in matural_transcript
    start_codon = matural_transcript.index(coding_sequence)
    # stop_codon means the last base position in matural transcript
    cds_len = len(coding_sequence)
    stop_codon = start_codon + cds_len
    cds_intervel = str(start_codon) + '-' + str(stop_codon)
    mt_len = len(matural_transcript)
    utr5 = matural_transcript[:start_codon]
    utr5_len = len(utr5)
    utr3 = matural_transcript[stop_codon:]
    utr3_len = len(utr3)
    for i in range(0, utr3_len-3):
        # start codon find 
        if utr3[i:i+3] == "ATG":
            for j in range(i+3, utr3_len, 3):
            # stop codon find
                if utr3[j:j+3] in stop_codon_list:
                    dORF = utr3[i:j+3]
                    out1 = [transcript_id, chrom, strand, cds_intervel,
                            stop_codon+i+1, stop_codon+j+3, 'dORF', len(dORF), dORF]
                    if not dORF_dict.get('dORF'):
                        dORF_dict['dORF'] = [out1]
                    else:
                        dORF_dict['dORF'].append(out1)
                    break 
    return dORF_dict



def get_mt_cds_seq(g, t, db, genome, style):
    """
    Extract transcipt and cds sequences and positions
    :param params:
    g (str): gene id
    t (str): transcript id
    db: db database from create
    genome: Fasta objective
    style (str): gtf/gff
    :return
    mt : str of transcript sequence
    cds: str of CDS sequence
    exons_list: List[Tuple(start, end)]
    cds_list: List[Tuple(start, end)]
    """
    # primary transcript (pt) 是基因组上的转录本的序列，
    # 有的会包括intron，所以提取的序列和matural transcript (mt) 不一致
    # pt = t.sequence(genome, use_strand=True)
    # matural transcript (mt)
    # 该步骤是获取序列: 转录本序列，CDS序列, five_prime_UTR, three_prime_UTR
    # 但是有的基因注释不正确，导致成熟转录本序列的两个来源不一致，一个是exon来源的，一个是five_prime_UTR和CDS来源的
    # 针对这种请情况，我们使用 five_prime_UTR + CDS + three_prime_UTR 拼接形成成熟的转录本序列用于uORF/dORF的注释
    # exon 提取的是转录本内成熟的mRNA的序列,即外显子按基因组坐标首尾相连的matural transcript
    mt, cds = '', ''
    exons_list = []
    mRNA = db[g][t]['mRNA']
    if db[g][t].get('exon'):
        for e in sorted(db[g][t]['exon'], key=lambda x : x.start):
            #不管是+还是-，一律从小到大排列，（有些GFF负链是从大到小排列的）
            exons_list.append((e.start, e.end))
            chrom, start, end, strand = e.chr, e.start, e.end, e.strand
            mt += str(genome[chrom][start-1 : end])
    mt = Seq(mt)
    if mRNA.strand == '-':
        mt = mt.reverse_complement()
    # CDS 提取的是编码区 ATG ... TAG
    cds_list = []
    if db[g][t].get('CDS'):
        for c in sorted(db[g][t]['CDS'], key=lambda x : x.start):
            #不管是+还是-，一律从小到大排列，（有些GFF负链是从大到小排列的）
            cds_list.append((c.start, c.end))
            chrom, start, end , strand = c.chr, c.start, c.end, c.strand
            cds += str(genome[c.chr][start-1:end])
    # CDS don't contains stop codon in gtf
    if db[g][t].get('CDS') and style == 'gtf':
        sc_seq = stop_codon_seq(db, g, t, genome[mRNA.chr])
        cds = add_stop_codon(cds, mRNA.strand, sc_seq)

    cds = Seq(cds)
    try:
        if len(cds) % 3 != 0:
            raise ValueError(f"The cds length of {t} ({len(cds)}) is not an integer multiple of 3. You need to check the annotation file.")
    except Exception as e:
        logging.warning(f"{e}")

    if mRNA.strand == '-':
        cds = cds.reverse_complement()
    mt = str(mt.upper())
    cds = str(cds.upper())

    return mt, cds, exons_list, cds_list


def dorf_filter_output(params:Tuple):
    """
    dorf filter by length threshold
    :param params:
    g (str): gene str
    t (str): transcript str
    db: db database from create
    genome: Fasta objective
    style (str): gtf/gff
    length (int): dORF length threshold
    output_format (str): output format (csv, gff, fasta)

    :return:
    dORF_seq_list: list object
    """
    g, symbol, t, db, genome, style, length, output_format = params
    index = 0
    # dORF_seq_list = deque()
    dORF_seq_list = list()
    # 可以多进程,但是多进程会导致和SQLite数据库的交互产生错误
    mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, style)
    mRNA = db[g][t]['mRNA']
    # objected-oriented
    if mt == cds: # not contain dORF
        return
    if len(cds) % 3 != 0: # CDS error
        return
    dORF_object = dORF(g, t, length, mRNA.chr, mRNA.strand, mt, cds)
    dORF_dict = dORF_object.dorf_parse()
    if len(dORF_dict) == 0: # CDS error
        return
    # loop for type1, type2 and type3
    for key in sorted(dORF_dict.keys()):
        # print(key,len(dORF_dict[key]))
        for it in dORF_dict[key]:
            # dORF length filter 
            if len(it[9]) <= length:
                # default: 6
                continue
            # csv output format 
            if output_format == 'fasta':
                dorf_id = ".".join(map(str,[it[1], it[5], it[6], it[7]]))
                # 'Gene ID', 'Transcript ID', 'Chrom', 'Strand',
                # 'CDS Interval', 'dORF Start', 'dORF End',
                # 'dORF Type', 'dORF Length', 'dORF'
                seq = it[9] # it[9] is str type
                desc='strand:%s dORF=%d-%d %s length=%d CDS=%s-%s'%(it[3],
                                                              it[5],
                                                              it[6],
                                                              it[7],
                                                              len(seq),
                                                              it[4].split('-')[0],
                                                              it[4].split('-')[1])
                dorfRecord = SeqRecord(Seq(seq), id=dorf_id.replace('transcript:',''), description=desc)
                dORF_seq_list.append(dorfRecord)

            elif output_format == 'gff':
                gff = GFF(exons_list, it, 'dORF')
                dorf_gff_line,features_gff_lines = gff.parse()
                dORF_seq_list.append(dorf_gff_line)
                for feature_line in features_gff_lines:
                    dORF_seq_list.append(feature_line)
            else:
                # csv output_format
                csv_row = it[:1] + [symbol] + it[1:]
                dORF_seq_list.append(dict(zip(_dORF_csv_header, csv_row)))
                index += 1
    return dORF_seq_list

def dorf_filter_output_faster(params:Tuple):
    """
    dorf filter by length threshold
    :param params:
    g (str): gene str
    t (str): transcript str
    db: db database from create
    genome: Fasta objective
    style (str): gtf/gff
    length (int): dORF length threshold
    output_format (str): output format (csv, gff, fasta)

    :return:
    dORF_seq_list: list object
    """
    g, symbol, t, db, genome, style, length, output_format = params
    index = 0
    # dORF_seq_list = deque()
    dORF_seq_list = list()
    # 可以多进程,但是多进程会导致和SQLite数据库的交互产生错误
    mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, style)
    mRNA = db[g][t]['mRNA']
    # objected-oriented
    if mt == cds: # not contain dORF
        return
    if len(cds) % 3 != 0: # CDS error
        return
    dORF_object = dORF(g, t, length, mRNA.chr, mRNA.strand, mt, cds)
    dORF_dict = dORF_object.dorf_parse()
    if len(dORF_dict) == 0: # CDS error
        return
    # loop for type1, type2 and type3
    for key in sorted(dORF_dict.keys()):
        # print(key,len(dORF_dict[key]))
        for it in dORF_dict[key]:
            # dORF length filter
            if len(it[9]) <= length:
                # default: 6
                continue
            # csv output format
            if output_format == 'fasta':
                dorf_id = ".".join(map(str,[it[1], it[5], it[6], it[7]]))
                # 'Gene ID', 'Transcript ID', 'Chrom', 'Strand',
                # 'CDS Interval', 'dORF Start', 'dORF End',
                # 'dORF Type', 'dORF Length', 'dORF'
                seq = it[9] # it[9] is str type
                desc='strand:%s dORF=%d-%d %s length=%d CDS=%s-%s'%(it[3],
                                                              it[5],
                                                              it[6],
                                                              it[7],
                                                              len(seq),
                                                              it[4].split('-')[0],
                                                              it[4].split('-')[1])
                dorfRecord = SeqRecord(Seq(seq), id=dorf_id.replace('transcript:',''), description=desc)
                dORF_seq_list.append(dorfRecord)

            elif output_format == 'gff':
                gff = GFF(exons_list, it, 'dORF')
                dorf_gff_line,features_gff_lines = gff.parse()
                dORF_seq_list.append(dorf_gff_line)
                for feature_line in features_gff_lines:
                    dORF_seq_list.append(feature_line)
            elif output_format == 'gtf':
                gff = GFF(exons_list, it, 'dORF')
                dorf_gff_line,features_gff_lines = gff.parse()
                # dORF
                dorf_gff_line_array = dorf_gff_line.split("\t")
                dorf_gff_line_array[8] = gff_attr_to_gtf(dorf_gff_line_array[8])
                dorf_gff_line_array[8] = dorf_gff_line_array[8].replace('ID', 'transcript_id')
                dorf_gff_line = "\t".join(dorf_gff_line_array)
                dORF_seq_list.append(dorf_gff_line)
                # CDS
                for feature_line in features_gff_lines:
                    feature_line_array = feature_line.split("\t")
                    feature_line_array[8] = gff_attr_to_gtf(feature_line_array[8])
                    feature_line_array[8] = feature_line_array[8].replace('Parent', 'transcript_id')
                    feature_line_array[8] = feature_line_array[8].replace('ID', 'dorf_id')
                    feature_line = "\t".join(feature_line_array)
                    dORF_seq_list.append(feature_line)
            else:
                # csv output_format
                dORF_seq_list.append(it[:1] + [symbol] + it[1:])
                index += 1
    return dORF_seq_list


def get_dorf(args):
    """
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout 
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="dORF.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    logging.info("Loading genome")
    genome = genome_dict(args.genome)

    # loop all transcripts in genome
    if not args.transcript:
        logging.info("dORF processing at whole-genome scale")
        param_iter = (
            (
                g, symbol, t, {g: {t: transcript}}, genome, args.style,
                args.length, args.output_format,
            )
            for g, symbol, t, transcript in iter_transcript_models_with_symbols(
                db, ('exon', 'CDS', 'stop_codon'),
                include_feature_attributes=False,
            )
        )
        # dORF_seq_list = deque()

        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_dORF_csv_header) + '\n')

        for para in tqdm(param_iter, ncols = 80, total=transcript_count(db, t2g), \
                         desc='dORF processing:', colour="GREEN"):
            dORF_seq_list = dorf_filter_output_faster(para)
            if dORF_seq_list == None:
                continue
            stream_parse_output(args, dORF_seq_list, out)
        # dORF_seq_list = [d for de in dORF_seq_list if de != None for d in de]
        # output file
        # parse_output(args, dORF_seq_list)
        out.close()
        logging.info("dORF finished at whole-genome scale")
    else:
        logging.info("dORF processing at transcript scale")
        t = args.transcript
        if not t2g.get(t):
            print('transcript id not exist in database, please check.')
            sys.exit(1)
        g = t2g[t]
        mRNA = db[g][t]['mRNA']
        para = (
            g, get_gene_symbol(db, g), t, db, genome, args.style,
            args.length, args.output_format,
        )
        dORF_seq_list = dorf_filter_output_faster(para)
        mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, args.style)
        dORF_object = dORF(g, t, args.length, mRNA.chr, mRNA.strand, mt, cds)

        dORF_dict = dORF_object.dorf_parse()
        if len(dORF_dict) == 0:
            print('The transcript does not contain dORF.')
            sys.exit(1)

        dorf_location_genome = [] # for schematic on genome # 3D list 
        for key in sorted(dORF_dict.keys()):
            for it in dORF_dict[key]:
                # dORF length filter 
                if len(it[9]) <= args.length:
                    # default: 6
                    continue
                # for schematic on genome
                gff = GFF(exons_list, it, 'dORF')
                ex_locations = gff.uorf_exons_location()
                dorf_location_genome.append(ex_locations)
        # figure the schametic 
        # without intron / -m 
        if args.schematic_without_intron:
            figure = visual_transcript(args.schematic_without_intron,
                                       args.transcript,
                                       dORF_object.transcript_location(),
                                       dORF_object.cds_location_transcript(),
                                       dORF_object.dorf_location_transcript(),
                                       'dORF')
            figure.draw()
        # with intron / -n 
        if args.schematic_with_intron:
            figure = visual_transcript_genome(mRNA.strand,
                                              args.schematic_with_intron,
                                              args.transcript,
                                              exons_list,
                                              cds_list,
                                              dorf_location_genome,# 3D list 
                                              'dORF')
            figure.draw()
        # file out
        parse_output(args, dORF_seq_list)
        logging.info("dORF finished at transcript scale")
