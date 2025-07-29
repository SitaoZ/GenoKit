# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from featurExtract.database.database import create_db, genome_dict
from featurExtract.utils.util import utr3_type, utr5_type, mRNA_type
from featurExtract.utils.util import gff_feature_dict, gtf_feature_dict
from featurExtract.utils.util import stop_codon_seq, add_stop_codon, parse_output
from featurExtract.utils.visualize import visual_transcript, visual_transcript_genome

# table header 
_uORF_csv_header = ['TranscriptID', 'Chrom', 'Strand', 
                    'CDS Interval', 'uORF Start', 'uORF End', 
                    'uORF Type', 'uORF Length', 'uORF']

class GFF(object):
    """A gff parse for line of uORF.
    
    Attributes:
        exons (2D list): [start,end] exon list in specifity transcript
        uorf_out (1D list): a list format including "TranscriptID,Chrom,
                         Strand,CDS Interval,uORF Start,uORF End,
                         uORF Type,uORF Length,uORF" from uORF.parse
    """
    def __init__(self, exons, uorf_out, feature_type):
        """ Initialize the usefull variable

        Args:
            exons (2D list): exon start and end list from a specifity transcript
            uorf_out (list): a list from uORF.parse 
            t_id (str) : transcript id
            chrom (str): chromosome id 
            strand (str) : transcript strand on genome
            start (str) : uORF start point on transcript
            end (str) : uORF end point on transcript
            uorf_type (str): uORF type 
            uorf_id (str): uORF id 
        """
        self.exons = exons
        self.uorf_out = uorf_out
        self.feature_type = feature_type
        self.t_id = self.uorf_out[0].replace('transcript:',self.feature_type+".")
        self.chrom = self.uorf_out[1] 
        self.strand = self.uorf_out[2]
        self.start = self.uorf_out[4]
        self.end = self.uorf_out[5]
        self.uorf_type = self.uorf_out[6]
        if self.feature_type == 'uORF':
            self.uorf_id = ".".join(map(str,[self.t_id, self.start, self.end, self.uorf_type]))
        else:
            self.uorf_id = ".".join(map(str,[self.t_id, self.start, self.end]))
    def parse(self):
        """ Parse the uORF into GFF line"""
        cds_interval = self.uorf_out[3]
        uORF_start, uORF_end, uORF_type = self.uorf_out[4:7]
        # build a class for uROF genome location 
        ugl = uORF_genome_location(self.exons, uORF_start, uORF_end, self.strand, cds_interval)
        genome_start, genome_end = ugl.parse()
        features = ugl.split_feature() # uORF exons blocks, genome location 
        features_gff_lines = self.feature_gff_format(features)
        #if strand == '-':
            # exchange location, because start less then end in gff 
        #    genome_start, genome_end = genome_end, genome_start
        # gff annotation line https://genome.ucsc.edu/FAQ/FAQformat#format3
        # col9 : id start end type 
        uorf_gff_line = "\t".join([self.chrom, "featurExtract", self.feature_type, \
                   str(genome_start), str(genome_end), ".", self.strand, ".", "ID="+self.uorf_id])
        
        return uorf_gff_line, features_gff_lines
        
    def uorf_exons_location(self):
        '''
        extract uorf exons blocks in genome, no reverse 
        '''
        cds_interval = self.uorf_out[3]
        uORF_start, uORF_end, uORF_type = self.uorf_out[4:7]
        # build a class for uROF genome location 
        ugl = uORF_genome_location(self.exons, uORF_start, uORF_end, self.strand, cds_interval)
        genome_start, genome_end = ugl.parse()
        features = ugl.split_feature()
        return features
    
    def parse_frame(self, length, frame):
        '''phase CDS: 0, 1, 2
           lenght: int,   previous feature length
           frame: (0,1,2),previous feature frame
        '''
        return (3- ((length-frame) % 3)) % 3
    
    def feature_gff_format(self, features):
        '''
        formate the split feature into gff line
        '''
        gff_lines = []
        frames = []
        lengths = []
        if self.strand == "+":
            for i,exon in enumerate(features):
                # http://mblab.wustl.edu/GTF22.html
                lengths.append(exon[1]-exon[0]+1)
                if i == 0:
                    frame = 0
                else:
                    frame = self.parse_frame(lengths[i-1], frames[i-1])
                frames.append(frame)
                gff_9_col = "ID="+self.uorf_id+":cds:"+str(i+1)+";Parent="+self.uorf_id
                gff_line = "\t".join([self.chrom, "featurExtract", "CDS", \
                   str(exon[0]), str(exon[1]), ".", self.strand, str(frame), gff_9_col])
                gff_lines.append(gff_line)
            return gff_lines
        else:
            # gff annotation line https://genome.ucsc.edu/FAQ/FAQformat#format3
            # features = features[::-1] # reverse, 位置顺序排列
            for i,exon in enumerate(features):
                lengths.append(exon[1]-exon[0]+1)
                if i == 0:
                    frame = 0
                else:
                    frame = self.parse_frame(lengths[i-1], frames[i-1])
                frames.append(frame)
            for i,exon in enumerate(features):
                cds_index = len(features) - i 
                gff_9_col = "ID="+self.uorf_id+":cds:"+str(cds_index)+";Parent="+self.uorf_id
                gff_line = "\t".join([self.chrom, "featurExtract", "CDS", \
                   str(exon[0]), str(exon[1]), ".", self.strand, str(frames[i]), gff_9_col])
                gff_lines.append(gff_line)
            return gff_lines[::-1]
        

class uORF_genome_location(object):
    """A parse for uORF location from transcript map to genome 
    designed for GFF output format 
    uORF的起始和终止只会位于某一外显子中,锁定外显子后根据相对位置判定真实坐标,
    其中uORF的起始等于它位于的外显子的起始位置+该位置离该exon起始的距离,
    uORF的终止等于它位于的外显子的起始位置+该位置离该exon起始的距离.
    """
    def __init__(self, exons, uorf_start, uorf_end, strand, cds_interval):
        self.exons = exons
        self.uorf_start = uorf_start
        self.uorf_end = uorf_end
        self.strand = strand 
        # strand affect the utr5 location, reverse -
        if self.strand == '-':
            self.exons = self.exons[::-1]
        self.cds_interval = cds_interval # for 5utr and 3utr
        if self.cds_interval:
            cdss = self.cds_interval.split('-')
            self.cds_start = cdss[0]
            self.cds_end = cdss[1]
    def parse(self):
        """uORF transcript location to genome location (start ,end)"""
        
        genome_start = self.location(self.uorf_start)
        genome_end = self.location(self.uorf_end)
        if self.strand == '-':
            # exchange location, because start less then end in gff 
            genome_start, genome_end = genome_end, genome_start 
        return genome_start, genome_end
        
    def location(self, position):
        """Used for uORF start and end mapping to genome location  
        
        Params:
            - position: the uORF start and end position in cDNA(mRNA)
        Return:
            - the uORF start or uORF end posiiton in genome 
        """ 
        exon_len_sum = 0
        genome_location = 0
        for i, exon in enumerate(self.exons):
            exon_len = exon[1] - exon[0] + 1 # 1 based; length should add 1
            prior_exon_len = exon_len_sum    # pay attention to the coordinate, assignment should be prior to sum
            exon_len_sum += exon_len 
            if position <= prior_exon_len:
                # uORF起始于第一个外显子
                if self.strand == '+':
                    genome_location = exon[0] + position - 1 # 1 based
                else:
                    genome_location = exon[1] - position + 1
                return genome_location
            # 长度位于前面累加和最新一个的中间
            elif prior_exon_len < position <= exon_len_sum:
                # UTR 存在内含子, uORF 起始于非第一个外显子
                # length 非第一个外显子上且距离所在外显子起始位置的距离
                length = position - prior_exon_len
                if self.strand == '+':
                    genome_location = exon[0] + length - 1 # 1 based 
                else:
                    genome_location = exon[1] - length + 1
                return genome_location
    def exon_index(self, position):
        """Position locate the exon in mRNA
        Params:
             - position: the index in cDNA(mRNA)
        Return:
             - the exon index containing position
        """
        exon_len_sum = 0
        for i, exon in enumerate(self.exons):
            exon_len = exon[1] - exon[0] + 1 # 1 based; length should add 1
            prior_exon_len = exon_len_sum    # pay attention to the coordinate, assignment should be prior to sum
            exon_len_sum += exon_len
            if prior_exon_len < position <= exon_len_sum:
                return i
    def split_feature(self):
        """Split the uORF into feature including exon and intron"""
        uorf_g_start, uorf_g_end = self.parse()
        # true position in genome
        if self.strand == '-':
            uorf_g_start, uorf_g_end = uorf_g_end, uorf_g_start
        u_start_exon_index, u_end_exon_index = self.exon_index(self.uorf_start), self.exon_index(self.uorf_end)
        # cds_g_start, cds_g_end = self.location((self.cds_start), self.location(self.cds_end)
        if u_start_exon_index == u_end_exon_index:
            # uORF start and end in same exon, no intron
            # type1, type2, type3 uORF
            if self.strand == '+':
                return [[uorf_g_start,uorf_g_end]]
            else:
                return [[uorf_g_end,uorf_g_start]]
        elif u_start_exon_index < u_end_exon_index:
            # have at least once intron, CDS feature collect
            # type1, type2, type3 uORF
            if self.strand == '+':
                s = [[uorf_g_start, self.exons[u_start_exon_index][1]]]
                e = [[self.exons[u_end_exon_index][0], uorf_g_end]]
            else:
                s = [[self.exons[u_start_exon_index][0], uorf_g_start]]
                e = [[uorf_g_end, self.exons[u_end_exon_index][1]]]
            
            inters = []
            for i,inter in enumerate(self.exons[u_start_exon_index+1:u_end_exon_index]):
                if len(inter) != 0:
                    inters.append(inter)
            features = s + inters + e # 从大到小
            return features
        else:
            sys.exit('exon index error')

class uORF(object):
    """uORF finder using object-oriented
    Attributes:
        transcript_id (str): transcript id
        chorm (str): chromosome id 
        strand (str): transcript strand 
        matural_transcript (Seq of Biopython): matural transcript sequence
        coding_sequence (Seq of Biopython): conding sequence of specifity transcript 
    """
    START_CODON = 'ATG'
    STOP_CODON_LIST = ['TAA','TAG','TGA']
    
    def __init__(self, transcript_id, chrom, strand, matural_transcript, coding_sequence):
        self.transcript_id = transcript_id # transcript id, str
        self.chrom = chrom                 # chromosome id ,str
        self.strand = strand               # transcript strand, str
        self.mt = matural_transcript       # a Seq type (Biopython) from mature transcript without intron
        self.cds = coding_sequence         # a Seq type (Biopython) from coding sequence start with ATG ; end with TAA, TGA, TAG
    def transcript_location(self):
        """Transcript start and end
        1D list regradless of strand  
        """    
        return [0,len(self.mt)]
    def cds_location_transcript(self):
        """Extract the cds start and end in transcript
        1D list regardless of strand 
        """
        start_codon_position = self.mt.index(self.cds)
        cds_start = start_codon_position + 1
        cds_end = start_codon_position + len(self.cds)
        return [cds_start, cds_end]
        
    def uorf_location_transcript(self):
        """Extract uorfs start and end in transcript
        2D list regardless strand 
        """ 
        d2 = []
        uORF_dict = self.uorf_parse()
        for uorf_type in uORF_dict.keys():
            for it in uORF_dict[uorf_type]:
                uorf_start, uorf_end = it[4], it[5]
                d2.append([uorf_start, uorf_end])
        return d2 
        
    def uorf_parse(self):
        """uorf finder main program, return a dict"""
        # uORF_dict = defaultdict(list)
        uORF_dict = defaultdict()
        # start_codon_position means the first base position in matural_transcript
        if self.cds not in self.mt:
            return uORF_dict
        start_codon_position = self.mt.index(self.cds)
        # stop_codon_position means the last base position in matural transcript
        stop_codon_position = start_codon_position + len(self.cds)
        # CDS on mt coordinate
        cds_interval = str(start_codon_position + 1) + '-' + str(stop_codon_position) # transcript coordinate 1 based 
        mt_len = len(self.mt)
        cds_len = len(self.cds)
        utr5 = self.mt[:start_codon_position]
        utr5_len = len(utr5)
        utr3 = self.mt[stop_codon_position:]
        utr3_len = len(utr3)
        
        # UTR5 range for ATG searching
        for i in range(0, utr5_len-3+1):
            # start codon find 
            if self.mt[i:i+3] == uORF.START_CODON:
                for j in range(i+3, stop_codon_position, 3):
                # stop codon find
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j < utr5_len:
                        # type1 uORF  upstream; not unique 
                        type1_uORF = self.mt[i:j+3]
                        out1 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type1', len(type1_uORF), type1_uORF]
                        if not uORF_dict.get('type1_uORF'):
                            uORF_dict['type1_uORF'] = deque([out1])
                        else:
                            uORF_dict['type1_uORF'].append(out1)
                        # only use the first stop codon, so break must
                        break
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j + 3 > utr5_len and j + 3 < stop_codon_position:
                        # type2 uORF across; the overlap region is triple or not; not unique
                        type2_uORF = self.mt[i:j+3]
                        out2 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type2', len(type2_uORF), type2_uORF]
                        if not uORF_dict.get('type2_uORF'):
                            uORF_dict['type2_uORF'] = deque([out2])
                        else:
                            uORF_dict['type2_uORF'].append(out2)
                        break
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j + 3 > utr5_len and j+3 == stop_codon_position and (utr5_len - i)%3 == 0:
                        # type3 uORF; N extention 通读; not unique 
                        type3_uORF = self.mt[i:utr5_len+cds_len]
                        out3 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type3', len(type3_uORF), type3_uORF]
                        if not uORF_dict.get('type3_uORF'):
                            uORF_dict['type3_uORF'] = deque([out3])
                        else:
                            uORF_dict['type3_uORF'].append(out3)
                        break
        return uORF_dict
            


def uorf_procedure_oriented(transcript_id, chrom, strand, matural_transcript, coding_sequence):
    """uorf finder using procedure-oriented
    parameters:
     transcript_id:      transcript id
     matural_transcript: a Seq type (Biopython) from mature transcript without intron
     coding_sequence:    a Seq type (Biopython) from coding sequence start with ATG , 
                         end with TAA, TGA, TAG
    return:
     upper stream open reading frame
    """
    uORF_dict = defaultdict(list)
    stop_codon_list = ['TAA','TAG','TGA']
    # start_codon means the first base position in matural_transcript
    start_codon = matural_transcript.index(coding_sequence)
    # print('xxxx',start_codon)
    # stop_codon means the last base position in matural transcript
    cds_len = len(coding_sequence)
    stop_codon = start_codon + cds_len
    # print('yyyy',cds_len)
    cds_interval = str(start_codon+1) + '-' + str(stop_codon) # transcript coordinate 1 based 
    mt_len = len(matural_transcript)
    utr5 = matural_transcript[:start_codon]
    utr5_len = len(utr5)
    utr3 = matural_transcript[stop_codon:]
    # UTR5 range 
    for i in range(0, utr5_len-3+1):
        # start codon find 
        if matural_transcript[i:i+3] == "ATG":
            for j in range(i+3,stop_codon,3):
            # stop codon find
                if matural_transcript[j:j+3] in stop_codon_list and j < utr5_len:
                    # type1 uORF  upstream; not unique 
                    type1_uORF = matural_transcript[i:j+3]
                    out1 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type1', len(type1_uORF), type1_uORF]
                    if not uORF_dict.get('type1_uORF'):
                        uORF_dict['type1_uORF'] = [out1]
                    else:
                        uORF_dict['type1_uORF'].append(out1)
                    break 
                if matural_transcript[j:j+3] in stop_codon_list and j + 3 > utr5_len and j + 3 < stop_codon:
                    # type2 uORF across; the overlap region is triple or not; not unique
                    type2_uORF = matural_transcript[i:j+3]
                    out2 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type2', len(type2_uORF), type2_uORF]
                    if not uORF_dict.get('type2_uORF'):
                        uORF_dict['type2_uORF'] = [out2]
                    else:
                        uORF_dict['type2_uORF'].append(out2)
                    break 
                if matural_transcript[j:j+3] in stop_codon_list and j + 3 > utr5_len and j+3 == stop_codon and (utr5_len - i)%3 == 0:
                    # N extention 通读; not unique 
                    type3_uORF = matural_transcript[i:utr5_len+cds_len]
                    out3 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type3', len(type3_uORF), type3_uORF]
                    if not uORF_dict.get('type3_uORF'):
                        uORF_dict['type3_uORF'] = [out3]
                    else:
                        uORF_dict['type3_uORF'].append(out3)
                    break
    return uORF_dict


def get_mt_cds_seq(g, t, db, genome, style):
    """
    parameters:
        g: gene
        t: transcript
        genome: geneome fasta
        style: GTF/GFF
    """
    # primary transcript (pt) 是基因组上的转录本的序列，
    # 有的会包括intron，所以提取的序列和matural transcript (mt) 不一致
    # pt = t.sequence(genome, use_strand=True)
    # matural transcript (mt)
    # 该步骤是获取序列: 转录本序列，CDS序列, five_prime_UTR, three_prime_UTR
    # 但是有的基因注释不正确，导致成熟转录本序列的两个来源不一致，一个是exon来源的，一个是five_prime_UTR和CDS来源的
    # 针对这种请情况，我们使用 five_prime_UTR + CDS + three_prime_UTR 拼接形成成熟的转录本序列用于uORF的注释
    # exon 提取的是转录本内成熟的mRNA的序列,即外显子按基因组坐标首尾相连的matural transcript
    mt, cds = '', ''
    exons_list = []
    mRNA = db[g][t]['mRNA']
    if db[g][t].get('exon'):
        for e in db[g][t]['exon']:
            exons_list.append((e.start, e.end))
            chrom, start, end , strand = e.chr, e.start, e.end, e.strand
            mt += str(genome[chrom][start-1:end])
    mt = Seq(mt)
    if mRNA.strand == '-':
        mt = mt.reverse_complement()
    # CDS 提取的是编码区 ATG ... TAG
    cds_list = []
    if db[g][t].get('CDS'):
        for c in db[g][t]['CDS']:
            chrom, start, end , strand = c.chr, c.start, c.end, c.strand
            cds_list.append((c.start, c.end))
            cds += str(genome[c.chr][start-1:end])
    #if style == 'GTF':
    #    sc_seq = stop_codon_seq(db, g, t, genome)
    #    cds = add_stop_codon(cds, t.strand, sc_seq)
    cds = Seq(cds)
    ## 需要确定exon和mRNA的位置是否一致，不一致以mRNA为准
    # mt from utr5/3 and cds
    if mRNA.strand == '-':
        cds = cds.reverse_complement()
    return mt, cds, exons_list, cds_list
    
def uorf_filter_output(params):
    g, t, db, genome, style, length, output_format = params
    index = 0 
    uORF_seq_list = deque()
    mt, cds, exons_list, cds_list = get_mt_cds_seq(g, t, db, genome, style)
    mRNA = db[g][t]['mRNA']
    # objected-oriented
    if mt == cds: # not contain uORF
        return 
    if len(cds) % 3 != 0: # CDS error
        return 
    uORF_ = uORF(t, mRNA.chr, mRNA.strand, mt, cds)
    uORF_dict = uORF_.uorf_parse()
    if len(uORF_dict) == 0: # CDS error
        return 
    # loop for type1, type2 and type3
    for key in sorted(uORF_dict.keys()):
        # print(key,len(uORF_dict[key]))
        for it in uORF_dict[key]:
            # uORF length filter 
            if len(it[8]) <= length:
                # default: 6
                continue
            # csv output format 
            if output_format == 'fasta':
                uorf_id = ".".join(map(str,[it[0], it[4], it[5], it[6]]))
                seq = it[8] # it[8] is s Seq type
                desc='strand:%s uORF=%d-%d %s length=%d CDS=%s-%s'%(it[2],
                                                              it[4],
                                                              it[5],
                                                              it[6],
                                                              len(seq),
                                                              it[3].split('-')[0],
                                                              it[3].split('-')[1])
                uorfRecord = SeqRecord(seq, id=uorf_id.replace('transcript:',''), description=desc)
                uORF_seq_list.append(uorfRecord)

            elif output_format == 'gff':
                gff = GFF(exons_list, it, 'uORF')
                uorf_gff_line, features_gff_lines = gff.parse()
                uORF_seq_list.append(uorf_gff_line)
                for feature_line in features_gff_lines:
                    uORF_seq_list.append(feature_line)
            else:
                # csv output_format
                uORF_seq_list.append(dict(zip(_uORF_csv_header, it)))
                index += 1
    return uORF_seq_list

def get_uorf(args):
    """ Get uORF main command line entry
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    """
    #if args.style == 'GFF':
    #    db, t2g = gff_feature_dict(args.database, args.style)
    #else:
    #    db, t2g = gtf_feature_dict(args.database, args.style)
    with open(args.database, 'rb') as f:
        db, t2g = cPickle.load(f)
    all_transcripts = [t for g in db for t in db[g] if t != 'gene']
    genome = genome_dict(args.genome)
    if not args.transcript:
        #param_list = [(t2g[t], t, db, genome, args.style, args.length, args.output_format) for t in all_transcripts]
        #uORF_seq_list = deque()
        #for para in tqdm(param_list, ncols = 80, total=len(param_list), desc='uORF processing:'):
        #    uORF_seq_list.append(uorf_filter_output(para))
        #uORF_seq_list = [d for de in uORF_seq_list if de != None for d in de]
        # output file 
        #parse_output(args, uORF_seq_list)
        # 2
        param_list = [(t, args.database, args.genome, args.style, args.length, args.output_format)
                      for t in db.features_of_type(mRNA_str, order_by='start')]
        with Pool(processes=args.process) as p:
            uORF_seq_list = list(tqdm(p.imap(uorf_filter_output, param_list),
                                 total=len(param_list),
                                 ncols = 80, desc='uORF processing:'))
        uORF_seq_list = [d for de in uORF_seq_list if de != None for d in de]
        # output file 
        parse_output(args, uORF_seq_list)

    else:
        g = t2g[args.transcript]
        mRNA = db[g][args.transcript]['mRNA']
        # specify the transcript id; only one transcript
        params = (g, args.transcript, db, genome, args.style, args.length, args.output_format)
        uORF_seq_list = uorf_filter_output(params)
        mt, cds, exons_list, cds_list = get_mt_cds_seq(g, args.transcript, db, genome, args.style)
        uORF_ = uORF(args.transcript, mRNA.chr, mRNA.strand, mt, cds)
        # figure the schametic 
        uORF_dict = uORF_.uorf_parse()
        uorf_location_genome = [] # for schematic on genome # 3D list 
        for key in sorted(uORF_dict.keys()):
            for it in uORF_dict[key]:
                # uORF length filter 
                if len(it[8]) <= args.length:
                    # default: 6
                    continue
                # for schematic on genome
                gff = GFF(exons_list, it, 'uORF')
                ex_locations = gff.uorf_exons_location()
                uorf_location_genome.append(ex_locations)
        # without intron 
        if args.schematic_without_intron == True:
            figure = visual_transcript(args.schematic_without_intron, 
                                       args.transcript, 
                                       uORF_.transcript_location(),
                                       uORF_.cds_location_transcript(),
                                       uORF_.uorf_location_transcript(),
                                       'uORF')
            figure.draw()
        # with intron
        if args.schematic_with_intron == True:
            figure = visual_transcript_genome(mRNA.strand, 
                                              args.schematic_with_intron, 
                                              args.transcript, 
                                              exons_list, 
                                              cds_list,
                                              uorf_location_genome,# 3D list 
                                              'uORF')
            figure.draw() 
        # file out
        parse_output(args, uORF_seq_list)
