# -*- coding: utf-8 -*-
'''
filename: build_db.py
date: 2024/11/22 15:28
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
import csv
import sys
import gzip
import os
# import polars as pl
import pandas as pd 
from tqdm import tqdm
from Bio import SeqIO
from io import StringIO
from typing import List, Dict, Any
from collections import defaultdict, namedtuple

# gffutils version
#def stop_codon_seq(db, transcript, genome):
#    s = ''
#    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
#        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
#    return s



def open_auto(filename, mode="rt"):
    """
    Open .gz automatically.
    Parameters
    ----------
    filename : str
        GFF file path
    mode : str
        Default: "rt"
    Returns
    -------
    file handle
    """
    filename = str(filename)

    if filename.endswith(".gz"):
        return gzip.open(filename, mode)

    return open(filename, mode)

def start_codon_seq(db, gene, transcript, genome):
    # 有的GTF文件中，转录本存在CDS和stop_codon, 但是没有start_codon属性
    # 会导致start_codon 为空字典{}
    seq = ''
    start_codon = db[gene][transcript]['start_codon']
    if len(start_codon) != 0:
        for s in start_codon:
            start, end = s.start, s.end
            seq += str(genome[start-1 : end])
        return seq
    else:
        return ''

def stop_codon_seq(db, gene, transcript, genome):
    # 有的GTF文件中，转录本存在CDS和start_codon, 但是没有stop_codon属性
    # 会导致stop_codon 为空字典{}
    seq = ''
    stop_codon = db[gene][transcript]['stop_codon']
    stop_codon_sort = sorted(stop_codon, key=lambda x: x.start)
    if len(stop_codon) != 0:
        for s in stop_codon_sort:
            start, end = s.start, s.end
            seq += str(genome[start-1 : end])
        return seq
    else:
        return ''

def add_stop_codon(seq, strand, stop_codon_seq):
    '''
    use for gtf seq add stop codon in CDS extracing  
    parameter:
     seq: sequence string 
     strand : ranscript strand 
     stop_codon_seq : stop codong sequence (3 bases)
    return:
     whole CDS
    '''
    if strand == '+':
        seq += stop_codon_seq
    elif strand == '-':
        seq = stop_codon_seq + seq # 负链 stop codon 位置靠前
    else:
        # no clear strand info
        seq += stop_codon_seq
    return seq

def add_start_codon(seq, strand, start_codon_seq):
    '''

    :param seq:
    :param strand:
    :param start_codon_seq:
    :return:
    '''
    if strand == '+':
        seq = start_codon_seq + seq
    elif strand == '-':
        seq = seq + start_codon_seq# 负链 start codon 位置靠后
    else:
        # no clear strand info
        seq = start_codon_seq + seq
    return seq

def mRNA_type(file_format: str):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return: 
     mRNA type string
    '''
    file_format = file_format.lower()
    if file_format == 'gff':
        return 'mRNA'
        
    elif file_format == 'gtf':
        return 'transcript'
    else:
        sys.stderr.write("parameter --style/-s should be assign \n")
        sys.exit(1)

def utr3_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return: 
     utr3 type string
    '''
    if file_format == 'gff':
        return 'three_prime_UTR'
    elif file_format == 'gtf':
        return 'three_prime_utr'
        # return '3UTR'
    else:
        sys.stderr.write("parameter --style/-s should be assign \n")
        sys.exit(1)

def utr5_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return:
     utr5 type string
    '''
    if file_format == 'gff':
        return 'five_prime_UTR'
    elif file_format == 'gtf':
        return 'five_prime_utr'
        # return '5UTR'
    else:
        sys.stderr.write("parameter --style/-s should be assign \n")
        sys.exit(1)

def seq_upper_lower(seq,start,end):
    '''
    parameter:
     seq: sequence 
     start: start codon position 
     end: stop codon position 
    return:
     seq with upper and lower bases
    '''
    utr5 = seq[:start].lower()
    coding = seq[start:end].upper()
    utr3 = seq[end:].lower()
    return utr5+coding+utr3

def gene_symbol(attribute, gene_id=None):
    """
    Return a gene symbol from a GFF3/GTF ninth-column attribute string.

    GFF3 uses ``Name=...`` and GTF uses ``gene_name "..."``.  When the
    symbol is absent, return the supplied gene ID (or the ID parsed from the
    attribute) so every output row has a useful, non-empty identifier.
    """
    attribute = str(attribute or "")
    fallback = str(gene_id or "")

    name_match = re.search(r'(?:^|;)\s*gene_name\s+"([^"]+)"', attribute)
    if name_match:
        return name_match.group(1)

    if re.search(r'(?:^|;)\s*(?:ID|Name)=', attribute):
        values = {}
        for field in attribute.split(";"):
            if "=" not in field:
                continue
            key, value = field.split("=", 1)
            values.setdefault(key.strip(), value.strip())
        return values.get("Name") or fallback or values.get("ID", "")

    id_match = re.search(r'(?:^|;)\s*gene_id\s+"([^"]+)"', attribute)
    return fallback or (id_match.group(1) if id_match else "")


def write_csv_csv(data_list: List[List[Any]],
                     output_path:str):
    """
    使用csv.writer写入（比DictWriter快15-30%）

    参数:
        example: 列表的列表 [[row1], [row2], ...]
        filepath: 输出文件路径
        headers: 表头列表
    """
    if not data_list:
        return

    with open(output_path, 'a+', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # 写入表头
        if isinstance(data_list[0], dict):
            headers = list(data_list[0].keys())
            writer.writerow(headers)
            data = [list(row.values()) for row in data_list]
        else:
            data = data_list

        # 批量写入
        writer.writerows(data)


def write_csv_polars(data: List[Dict[str, Any]],
                     filepath: str,
                     optimize: bool = True) -> None:
    """
    使用Polars写入CSV（性能最优）

    安装: pip install polars
    """
    # 1. 创建DataFrame（比pandas快5-10倍）
    df = pl.DataFrame(data)

    # 2. 写入CSV（支持多线程和优化）
    write_options = {
        'separator': ',',
        'include_header': True,
        'batch_size': 1024 * 1024,  # 1MB批处理
    }

    if optimize:
        write_options.update({
            'datetime_format': '%Y-%m-%d %H:%M:%S',
            'float_precision': 6,
        })

    df.write_csv(filepath, **write_options)



def parse_output(args, seq_list):
    """
    function for output
    parameters:
        args: a class for argparse
        seq_list: Promoter/gene/transcript/exon/intron/cds/utr/uORF/dORF/CDS seq_list
    return: none
    """
    if args.print:
        if args.output_format == 'csv':
            df_seq = pd.DataFrame.from_dict(seq_list)
            df_stream = StringIO()
            df_seq.to_csv(df_stream, sep=',', index=False)
            print(df_stream.getvalue(), file=sys.stdout, end="")
        elif args.output_format == 'fasta':
            SeqIO.write(seq_list, sys.stdout, "fasta")
        elif args.output_format == 'gff':
            for record in seq_list:
                print(record, file=sys.stdout)
        elif args.output_format == 'gtf':
            for record in seq_list:
                print(record, file=sys.stdout)

        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff, gtf}.')
    elif args.output != None:
        if args.output_format == 'csv':
            # 1.pandas
            df_seq = pd.DataFrame.from_dict(seq_list)
            df_seq.to_csv(args.output, sep=',', index=False)
            # 2.csv
            # write_csv_csv(seq_list, args.output)
            # 3.polars
            # write_csv_polars(seq_list, args.output)
        elif args.output_format == 'fasta':
            with open(args.output, 'w') as handle:
                SeqIO.write(seq_list, handle, "fasta")
        elif args.output_format == 'gff':
            with open(args.output, 'w') as handle:
                for record in seq_list:
                    handle.write(str(record)+'\n')
        elif args.output_format == 'gtf':
            with open(args.output, 'w') as handle:
                for record in seq_list:
                    handle.write(str(record)+'\n')
        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff, gtf}.')
    else:
        sys.exit('-o or -v should be specified')

def stream_parse_output(args, seq_list, handle):
    """
        function for output
        parameters:
            args: a class for argparse
            seq_list: Promoter/gene/transcript/exon/intron/cds/utr/uORF/dORF/CDS seq_list
        return: none
        """
    if args.print:
        if args.output_format == 'csv':
            df_seq = pd.DataFrame.from_dict(seq_list)
            df_stream = StringIO()
            df_seq.to_csv(df_stream, sep=',', index=False)
            print(df_stream.getvalue(), file=sys.stdout, end="")
        elif args.output_format == 'fasta':
            SeqIO.write(seq_list, sys.stdout, "fasta")
        elif args.output_format == 'gff':
            for record in seq_list:
                print(record, file=sys.stdout)
        elif args.output_format == 'gtf':
            for record in seq_list:
                print(record, file=sys.stdout)

        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff, gtf}.')
    elif args.output != None:
        if args.output_format == 'csv':
            for record in seq_list:
                line = ",".join(map(str,record))
                handle.write(line+'\n')
        elif args.output_format == 'fasta':
            SeqIO.write(seq_list, handle, "fasta")
        elif args.output_format == 'gff':
            for record in seq_list:
                handle.write(str(record) + '\n')
        elif args.output_format == 'gtf':
            for record in seq_list:
                handle.write(str(record) + '\n')
        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff, gtf}.')
    else:
        sys.exit('-o or -v should be specified')

def safe_count_lines(filename):
    """带异常处理的行数统计"""
    try:
        # gz 压缩文件的处理
        is_gzip = filename.lower().endswith("gz")
        opener = gzip.open if is_gzip else open
        with opener(filename, 'rb') as f:
            line_count = 0
            file_size = 0
            buf_size = 1024 * 1024
            
            while True:
                buf = f.read(buf_size)
                if not buf:
                    break
                file_size += len(buf)
                line_count += buf.count(b'\n')
            
            if file_size == 0:
                return 0
            if buf and buf[-1] != 10:  # 10是\n的ASCII码
                line_count += 1
            
            return line_count
    except FileNotFoundError:
        print(f"错误：文件 {filename} 不存在")
        return -1
    except Exception as e:
        print(f"读取文件时发生错误：{str(e)}")
        return -1

def gtf_gene(db, g):
    '''
    obtain gene info from its transcripts in gtf
    :param db:
    :param g:
    :return:
    '''
    transcripts = [
        transcript for key, transcript in db[g].items()
        if key not in {'gene', 'symbol'}
    ]
    chrom = transcripts[0]['mRNA'].chr
    source = transcripts[0]['mRNA'].source
    type = 'gene'
    start = min(transcript['mRNA'].start for transcript in transcripts)
    end = max(transcript['mRNA'].end for transcript in transcripts)
    score = ''
    strand = transcripts[0]['mRNA'].strand
    phase = '.'
    symbol = db[g].get('symbol') or g
    attribute = f'gene_id "{g}";'
    if symbol != g:
        attribute += f' gene_name "{symbol}";'
    array = [chrom, source, type, start, end, score, strand, phase, attribute]
    return record._make(array)

_RecordBase = namedtuple(
    'record',
    ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute'],
)


class record(_RecordBase):
    """Nine-column feature record with a derived gene-symbol property."""

    __slots__ = ()

    @property
    def symbol(self):
        return gene_symbol(self.attribute)

_GTF_GENE_ID = 'gene_id "'
_GTF_TRANSCRIPT_ID = 'transcript_id "'


def _quoted_gtf_attribute(attributes, prefix):
    """Return a standard quoted GTF attribute without regex overhead."""
    start = attributes.find(prefix)
    if start < 0:
        raise ValueError(f'Missing {prefix[:-1]} attribute in GTF record')
    start += len(prefix)
    end = attributes.find('"', start)
    if end < 0:
        raise ValueError(f'Unterminated {prefix[:-1]} attribute in GTF record')
    return attributes[start:end]

def gtf_feature_dict(file_path, file_type):
    """
    parse gtf genome feature file
      param: 
          file_path: gff/gtf path
          file_type: gff or gtf
      return:
          feature dict
          transcript2gene dict
    """
    five_prime = utr5_type(file_type)
    three_prime = utr3_type(file_type)
    g_dict = dict()
    t2g = dict()
    all_gene = set()
    # Count progress in bytes. This avoids a complete line-count pre-pass over
    # multi-gigabyte GTF files before parsing can begin.
    file_size = os.path.getsize(file_path)
    with open_auto(file_path, 'rb') as f, tqdm(
        total=file_size,
        unit='B',
        unit_scale=True,
        mininterval=0.5,
        ncols=80,
        desc='Reading annotation',
        colour='GREEN',
    ) as progress:
        for raw_line in f:
            progress.update(len(raw_line))
            if raw_line.startswith(b'#') or not raw_line.strip():
                continue
            line = raw_line.decode('utf-8').strip()
            array = line.split('\t', 8)
            f_type, array[3], array[4] = array[2], int(array[3]), int(array[4])
            # ensemble down gtf have gene feature
            if 'gene' in f_type:
                # some gtf do not has gene info
                # wget https://ftp.ensembl.org/pub/release-115/gtf/danio_rerio/Danio_rerio.GRCz11.115.gtf.gz
                # gtf 默认是外显子从小到大的顺序排 (hg38.ensGene.gtf)，但是ensemble的GTF负链的基因是从大道小排序(Danio_rerio.GRCz11.115.gtf.gz)
                # 要需要统一从大到小排序
                gene = _quoted_gtf_attribute(array[8], _GTF_GENE_ID)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene]['gene'] = record._make(array)     # record
                else:
                    if not g_dict[gene].get('gene'):
                        g_dict[gene]['gene'] = record._make(array) # record
                symbol = gene_symbol(array[8], gene)
                if not g_dict[gene].get('symbol') or g_dict[gene]['symbol'] == gene:
                    g_dict[gene]['symbol'] = symbol
                continue
            if f_type == 'miRNA':
                continue
            ts = _quoted_gtf_attribute(array[8], _GTF_TRANSCRIPT_ID)
            gene = _quoted_gtf_attribute(array[8], _GTF_GENE_ID)
            if 'RNA' in f_type or 'transcript' in f_type:
                t2g[ts] = gene
                all_gene.add(gene)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene][ts] = defaultdict(dict)
                else:
                    g_dict[gene][ts] = defaultdict(dict)
                symbol = gene_symbol(array[8], gene)
                if not g_dict[gene].get('symbol') or g_dict[gene]['symbol'] == gene:
                    g_dict[gene]['symbol'] = symbol
                g_dict[gene][ts]["mRNA"] = record._make(array) # record
                g_dict[gene][ts]["CDS"] = []
                g_dict[gene][ts]["exon"] = []
                g_dict[gene][ts][five_prime] = []
                g_dict[gene][ts][three_prime] = []
                g_dict[gene][ts]["start_codon"] = []
                g_dict[gene][ts]["stop_codon"] = []
            elif f_type == 'CDS':
                if gene not in all_gene:
                    continue
                g_dict[gene][ts]["CDS"].append(record._make(array))
            elif 'exon' in f_type:
                if gene not in all_gene:
                    continue
                g_dict[gene][ts]["exon"].append(record._make(array))
            elif f_type == five_prime:
                if gene not in all_gene:
                    continue
                g_dict[gene][ts][five_prime].append(record._make(array))
            elif f_type == three_prime:
                if gene not in all_gene:
                    continue
                g_dict[gene][ts][three_prime].append(record._make(array))
            elif f_type == 'start_codon':
                # 有些基因的start codon分布在不同的exon上面
                g_dict[gene][ts]['start_codon'].append(record._make(array))
            elif f_type == 'stop_codon':
                # 有些基因的stop codon分布在不同的exon上面
                g_dict[gene][ts]['stop_codon'].append(record._make(array))
            else:
                pass
    # check
    # g_dict format
    # g_dict: {'gene_id':
    #                    {'gene':record},
    #                    {'t1':{'mRNA':record, 'exon':[], 'CDS':[], 'five':[] , 'three':[], 'start_codon':[], 'stop_codon':[] }},
    #                    {'t2':{'mRNA':record, 'exon':[], 'CDS':[], 'five':[] , 'three':[], 'start_codon':[], 'stop_codon':[] }}
    #          }
    with tqdm(
        total=len(g_dict),
        unit='gene',
        mininterval=0.5,
        ncols=80,
        desc='Finalizing database',
        colour='GREEN',
    ) as progress:
        for g in g_dict:
            # add gene feature in case of the gtf containing gene info
            if not g_dict[g].get('gene'):
                g_dict[g]['gene'] = gtf_gene(g_dict, g)
            if not g_dict[g].get('symbol'):
                g_dict[g]['symbol'] = gene_symbol(g_dict[g]['gene'].attribute, g)
            for t in g_dict[g]:
                if t == 'gene':
                    continue
                if t == 'symbol':
                    continue
                if not g_dict[g][t].get('mRNA'):
                    g_dict[g][t]['mRNA'] = g_dict[g]['gene']
                # sort
                for feature in g_dict[g][t].values():
                    if not isinstance(feature, list) or len(feature) < 2:
                        continue
                    feature.sort(key=lambda item: item.start)
            progress.update(1)
    # print('have gene', g_dict['ENSG00000223972'])
    return g_dict, t2g


def gff_feature_dict(file_path, file_type):
    """
    parse gff genome feature file
      param: 
          file_path: gff/gtf path
          file_type: gff or gtf
      return:
          feature dict
          transcript2gene dict
    """
    t_pat = re.compile(r'ID=(\S+?);')
    g_pat = re.compile(r'Parent=(\S+?);')
    five_prime = utr5_type(file_type)
    three_prime = utr3_type(file_type)
    g_dict = dict()
    t2g = dict()
    all_gene = set()
    # Use byte progress to avoid a complete line-count pre-pass over large GFF
    # files.  This mirrors the GTF reader and starts parsing immediately.
    file_size = os.path.getsize(file_path)
    with open_auto(file_path, 'rb') as f, tqdm(
        total=file_size,
        unit='B',
        unit_scale=True,
        mininterval=0.5,
        ncols=80,
        desc='Reading annotation',
        colour='GREEN',
    ) as progress:
        for raw_line in f:
            progress.update(len(raw_line))
            if raw_line.startswith(b'#') or not raw_line.strip():
                continue
            line = raw_line.decode('utf-8').strip()
            if not line.endswith(";"):
                line = line + ';'
            array = line.split('\t')
            if len(array) < 5:
                print(array)
            array[3], array[4] = int(array[3]), int(array[4])
            f_type = array[2]
            if 'gene' in f_type:
                # gff has gene info
                gene = t_pat.search(array[8]).group(1)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene]['gene'] = record._make(array)
                else:
                    if not g_dict[gene].get('gene'):
                        g_dict[gene]['gene'] = record._make(array)
                symbol = gene_symbol(array[8], gene)
                if not g_dict[gene].get('symbol') or g_dict[gene]['symbol'] == gene:
                    g_dict[gene]['symbol'] = symbol
            if f_type == 'miRNA':
                continue
            if 'RNA' in f_type or 'transcript' in f_type:
                if 'gene' in f_type:
                    continue
                ts = t_pat.search(array[8]).group(1)
                gene = g_pat.search(array[8]).group(1)
                t2g[ts] = gene
                all_gene.add(gene)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene][ts] = defaultdict(dict)
                else:
                    g_dict[gene][ts] = defaultdict(dict)
                g_dict[gene][ts]["mRNA"] = record._make(array)
                g_dict[gene][ts]["CDS"] = []
                g_dict[gene][ts]["exon"] = []
                g_dict[gene][ts][five_prime] = []
                g_dict[gene][ts][three_prime] = []
                g_dict[gene][ts]["start_codon"] = []
                g_dict[gene][ts]["stop_codon"] = []
            elif f_type == 'CDS':
                tss = g_pat.search(array[8]).group(1)
                for ts in tss.split(','):
                    if g_dict.get(ts):
                        continue
                    if not t2g.get(ts):
                        continue
                    gene = t2g[ts]
                    if gene not in all_gene:
                        continue
                    g_dict[gene][ts]["CDS"].append(record._make(array))
            elif f_type == 'exon':
                tss = g_pat.search(array[8]).group(1)
                for ts in tss.split(','):
                    if g_dict.get(ts):
                        continue
                    if not t2g.get(ts):
                        continue
                    gene = t2g[ts]
                    if gene not in all_gene:
                        continue
                    if not g_dict[gene][ts].get('exon'):
                        g_dict[gene][ts]["exon"] = []
                    g_dict[gene][ts]["exon"].append(record._make(array))
            elif f_type == five_prime:
                tss = g_pat.search(array[8]).group(1)
                for ts in tss.split(','):
                    if not t2g.get(ts):
                        continue
                    gene = t2g[ts]
                    if gene not in all_gene:
                        continue
                    g_dict[gene][ts][five_prime].append(record._make(array))
            elif f_type == three_prime:
                tss = g_pat.search(array[8]).group(1)
                for ts in tss.split(','):
                    if not t2g.get(ts):
                        continue
                    gene = t2g[ts]
                    if gene not in all_gene:
                        continue
                    g_dict[gene][ts][three_prime].append(record._make(array))
            elif f_type == 'start_codon':
                # 有些基因的start codon分布在不同的exon上面
                g_dict[gene][ts]['start_codon'].append(record._make(array))
            elif f_type == 'stop_codon':
                # not a list
                g_dict[gene][ts]['stop_codon'].append(record._make(array))
            else:
                pass
    # Check and finalize in a visible second phase, matching GTF behavior.
    with tqdm(
        total=len(g_dict),
        unit='gene',
        mininterval=0.5,
        ncols=80,
        desc='Finalizing database',
        colour='GREEN',
    ) as progress:
        for g in g_dict:
            if not g_dict[g].get('symbol'):
                g_dict[g]['symbol'] = gene_symbol(g_dict[g]['gene'].attribute, g)
            for t in g_dict[g]:
                if t in {'gene', 'symbol'}:
                    continue
                if not g_dict[g][t].get('mRNA'):
                    g_dict[g][t]['mRNA'] = g_dict[g]['gene']
                # Most feature lists contain zero or one item.  Sorting only
                # longer lists avoids millions of temporary list allocations.
                for feature in g_dict[g][t].values():
                    if isinstance(feature, list) and len(feature) > 1:
                        feature.sort(key=lambda item: item.start)
            progress.update(1)

    return g_dict, t2g
