import re
import sys
import gzip 
import pandas as pd 
from tqdm import tqdm
from Bio import SeqIO
from io import StringIO
from collections import defaultdict, namedtuple

# gffutils version
#def stop_codon_seq(db, transcript, genome):
#    s = ''
#    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
#        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
#    return s


def stop_codon_seq(db, gene, transcript, genome):
    # 有的GFT文件中，转录本存在CDS和start_codon, 但是没有stop_codon属性
    # 会导致stop_codon 为空字典{}
    stop_codon = db[gene][transcript]['stop_codon']
    if len(stop_codon) != 0:
        start, end = stop_codon.start, stop_codon.end
        return str(genome[start-1 : end])
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

def mRNA_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return: 
     mRNA type string
    '''
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
        # return 'three_prime_utr'
        return 'three_prime_UTR'
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
        # return 'five_prime_utr'
        return 'five_prime_UTR'
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
            cds_seq = pd.DataFrame.from_dict(seq_list)
            cds_stream = StringIO()
            cds_seq.to_csv(cds_stream, sep=',', index=False)
            print(cds_stream.getvalue(), file=sys.stdout, end="")
        elif args.output_format == 'fasta':
            SeqIO.write(seq_list, sys.stdout, "fasta")
        elif args.output_format == 'gff':
            for record in seq_list:
                print(record, file=sys.stdout)
        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff}.')
    elif args.output != None:
        if args.output_format == 'csv':
            cds_seq = pd.DataFrame.from_dict(seq_list)
            cds_seq.to_csv(args.output, sep=',', index=False)
        elif args.output_format == 'fasta':
            with open(args.output, 'w') as handle:
                SeqIO.write(seq_list, handle, "fasta")
        elif args.output_format == 'gff':
            with open(args.output, 'w') as handle:
                for record in seq_list:
                    handle.write(str(record)+'\n')
        else:
            sys.exit('output(-o) format not be specified {csv, fasta, gff}.')
    else:
        sys.exit('-o or -v should be specified')


def safe_count_lines(filename):
    """带异常处理的行数统计"""
    try:
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


record = namedtuple('record', ['chr','source', 'type','start', 'end', 'score', 'strand', 'phase', 'attribute'])

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
    t_pat = re.compile(r'transcript_id "(\S+?)";') if file_type == 'gtf' else re.compile(r'ID=(\S+?);')
    g_pat = re.compile(r'gene_id "(\S+?)";') if file_type == 'gtf' else re.compile(r'Parent=(\S+?);')
    five_prime = utr5_type(file_type)
    three_prime = utr3_type(file_type)
    g_dict = dict()
    t2g = dict()
    all_gene = set()
    gene_info = dict()
    file_line_count = safe_count_lines(file_path)
    with open(file_path, 'r') as f:
        # for line in f:
        for line in tqdm(f, ncols = 80, total=file_line_count, desc='Database processing:'):
            if line.startswith('#'):
                continue
            line = line.strip()
            array = line.split('\t')
            array[3], array[4] = int(array[3]), int(array[4])
            f_type = array[2]
            if 'gene' in f_type:
                gene = g_pat.search(array[8]).group(1)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene]['gene'] = record._make(array)
                else:
                    if not g_dict[gene].get('gene'):
                        g_dict[gene]['gene'] = record._make(array)
                continue
            if f_type == 'miRNA':
                continue
            ts = t_pat.search(array[8]).group(1)
            gene = g_pat.search(array[8]).group(1)
            if 'RNA' in f_type or 'transcript' in f_type:
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
            elif f_type == 'stop_codon':
                g_dict[gene][ts]['stop_codon'] = record._make(array)
            else:
                pass
    # check
    for g in g_dict:
        for t in g_dict[g]:
            if t == 'gene':
                continue
            if not g_dict[g][t].get('mRNA'):
                g_dict[g][t]['mRNA'] = g_dict[g]['gene']
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
    file_line_count = safe_count_lines(file_path)
    with open(file_path, 'r') as f:
        # for line in f:
        for line in tqdm(f, ncols = 80, total=file_line_count, desc='Database processing:'):
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line.endswith(";"):
                line = line + ';'
            array = line.split('\t')
            array[3], array[4] = int(array[3]), int(array[4])
            f_type = array[2]
            if 'gene' in f_type:
                gene = t_pat.search(array[8]).group(1)
                if not g_dict.get(gene):
                    g_dict[gene] = defaultdict(dict)
                    g_dict[gene]['gene'] = record._make(array)
                else:
                    if not g_dict[gene].get('gene'):
                        g_dict[gene]['gene'] = record._make(array)
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
            elif f_type == 'stop_codon':
                g_dict[gene][ts]['stop_codon'] = record._make(array)
            else:
                pass
    # check
    for g in g_dict:
        for t in g_dict[g]:
            if t == 'gene':
                continue
            if not g_dict[g][t].get('mRNA'):
                g_dict[g][t]['mRNA'] = g_dict[g]['gene']
    # sort by position 
    # gtf_list.sort(key=lambda x: (x.seqname, x.attribute['gene_id'], x.attribute['transcript_id'], x.start, x.end))
    return g_dict, t2g

