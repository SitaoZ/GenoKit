import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, MeltingTemp
from Bio.Alphabet import generic_dna
from Bio.Emboss import PrimerSearch
import subprocess
import tempfile

class PrimerDesigner:
    def __init__(self, genome_file, annotation_file, output_dir="primers"):
        """
        初始化引物设计器
        
        参数:
            genome_file (str): 基因组序列文件路径(FASTA格式)
            annotation_file (str): 基因组注释文件路径(GFF/GTF格式)
            output_dir (str): 输出目录路径
        """
        self.genome_file = genome_file
        self.annotation_file = annotation_file
        self.output_dir = output_dir
        self.genome = None
        self.annotations = []
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 加载基因组和注释
        self._load_genome()
        self._load_annotations()
    
    def _load_genome(self):
        """加载基因组序列"""
        self.genome = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
    
    def _load_annotations(self):
        """加载基因组注释(GFF/GTF格式)"""
        with open(self.annotation_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                # 提取注释信息
                chrom = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                
                # 解析属性字段
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                    elif ' ' in attr:
                        key, value = attr.split(' ', 1)
                        attr_dict[key] = value.strip('"')
                
                # 只保留基因和转录本信息
                if feature_type in ['gene', 'mRNA', 'transcript']:
                    self.annotations.append({
                        'chrom': chrom,
                        'type': feature_type,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'attributes': attr_dict
                    })
    
    def get_gene_sequence(self, gene_id, upstream=500, downstream=500):
        """
        获取基因序列及其上下游区域
        
        参数:
            gene_id (str): 基因ID
            upstream (int): 上游区域长度(bp)
            downstream (int): 下游区域长度(bp)
            
        返回:
            Seq: 基因序列(包含上下游)
            dict: 基因注释信息
        """
        gene_info = None
        for ann in self.annotations:
            if ann['attributes'].get('ID', '').startswith(gene_id) or \
               ann['attributes'].get('gene_id', '').startswith(gene_id) or \
               ann['attributes'].get('Name', '').startswith(gene_id):
                gene_info = ann
                break
        
        if not gene_info:
            raise ValueError(f"Gene {gene_id} not found in annotations")
        
        chrom_seq = self.genome[gene_info['chrom']].seq
        
        # 计算序列区域
        start = max(1, gene_info['start'] - upstream)
        end = min(len(chrom_seq), gene_info['end'] + downstream)
        
        gene_seq = chrom_seq[start-1:end]
        
        # 如果是负链，取反向互补
        if gene_info['strand'] == '-':
            gene_seq = gene_seq.reverse_complement()
        
        return gene_seq, gene_info
    
    def design_primers(self, gene_id, primer_length=20, tm_min=55, tm_max=65, 
                       gc_min=40, gc_max=60, product_size=(100, 500), 
                       max_self_complementarity=5, max_3_self_complementarity=3):
        """
        设计基因特异性引物
        
        参数:
            gene_id (str): 基因ID
            primer_length (int): 引物长度
            tm_min (float): 最小熔解温度(°C)
            tm_max (float): 最大熔解温度(°C)
            gc_min (float): 最小GC含量(%)
            gc_max (float): 最大GC含量(%)
            product_size (tuple): 产物大小范围(bp)
            max_self_complementarity (int): 最大自身互补性
            max_3_self_complementarity (int): 3'端最大自身互补性
            
        返回:
            dict: 包含引物对和验证结果的字典
        """
        # 获取基因序列
        gene_seq, gene_info = self.get_gene_sequence(gene_id)
        gene_seq_str = str(gene_seq)
        
        # 使用EMBOSS的eprimer3设计引物
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta') as tmp_seq:
            # 写入临时序列文件
            tmp_seq.write(f">target_gene_{gene_id}\n{gene_seq_str}\n")
            tmp_seq.flush()
            
            # 准备eprimer3参数
            output_file = os.path.join(self.output_dir, f"{gene_id}_primers.out")
            args = [
                'eprimer3',
                '-sequence', tmp_seq.name,
                '-outfile', output_file,
                '-primerlength', str(primer_length),
                '-mintm', str(tm_min),
                '-maxtm', str(tm_max),
                '-maxgc', str(gc_max),
                '-mingc', str(gc_min),
                '-productosize', f"{product_size[0]}-{product_size[1]}",
                '-maxselfcomplementarity', str(max_self_complementarity),
                '-max3selfcomplementarity', str(max_3_self_complementarity),
                '-numreturn', '5'
            ]
            
            # 运行eprimer3
            subprocess.run(args, check=True)
            
            # 解析输出结果
            primers = self._parse_eprimer3_output(output_file)
        
        # 验证引物特异性
        if primers:
            for primer_pair in primers:
                self._check_specificity(primer_pair, gene_id)
        
        return {
            'gene_id': gene_id,
            'gene_info': gene_info,
            'primers': primers
        }
    
    def _parse_eprimer3_output(self, output_file):
        """解析eprimer3输出文件"""
        primers = []
        
        with open(output_file, 'r') as f:
            lines = f.readlines()
            
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith('PRIMER PAIR'):
                    pair_info = {}
                    
                    # 提取引物对信息
                    while i < len(lines) and not lines[i].startswith('='):
                        line = lines[i].strip()
                        
                        if line.startswith('PRIMER PAIR'):
                            pair_num = line.split()[2].strip(':')
                            pair_info['pair_num'] = pair_num
                        
                        elif line.startswith('PRODUCT SIZE:'):
                            pair_info['product_size'] = int(line.split(':')[1].strip())
                        
                        elif line.startswith('FORWARD PRIMER'):
                            parts = line.split()
                            pair_info['forward'] = {
                                'seq': parts[3],
                                'start': int(parts[5]),
                                'length': int(parts[7]),
                                'tm': float(parts[9]),
                                'gc': float(parts[11])
                            }
                        
                        elif line.startswith('REVERSE PRIMER'):
                            parts = line.split()
                            pair_info['reverse'] = {
                                'seq': parts[3],
                                'start': int(parts[5]),
                                'length': int(parts[7]),
                                'tm': float(parts[9]),
                                'gc': float(parts[11])
                            }
                        
                        i += 1
                    
                    primers.append(pair_info)
                else:
                    i += 1
        
        return primers
    
    def _check_specificity(self, primer_pair, target_gene_id):
        """
        检查引物特异性，确保不会扩增其他基因
        
        参数:
            primer_pair (dict): 引物对信息
            target_gene_id (str): 目标基因ID
        """
        # 创建临时文件用于引物搜索
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta') as tmp_primer:
            # 写入引物序列
            tmp_primer.write(f">forward_primer\n{primer_pair['forward']['seq']}\n")
            tmp_primer.write(f">reverse_primer\n{primer_pair['reverse']['seq']}\n")
            tmp_primer.flush()
            
            # 使用EMBOSS的primersearch检查特异性
            output_file = os.path.join(self.output_dir, f"{target_gene_id}_specificity.out")
            args = [
                'primersearch',
                '-seqall', self.genome_file,
                '-infile', tmp_primer.name,
                '-mismatchpercent', '20',  # 允许20%的错配
                '-outfile', output_file
            ]
            
            subprocess.run(args, check=True)
            
            # 解析结果
            self._parse_primersearch_output(output_file, primer_pair, target_gene_id)
    
    def _parse_primersearch_output(self, output_file, primer_pair, target_gene_id):
        """解析primersearch输出文件，检查特异性"""
        with open(output_file, 'r') as f:
            content = f.read()
            
            # 检查是否有多于一个的扩增产物
            if 'Amplimer count' in content:
                amplimer_counts = []
                for line in content.split('\n'):
                    if line.startswith('Amplimer count'):
                        count = int(line.split(':')[1].strip())
                        amplimer_counts.append(count)
                
                # 如果总扩增数大于1，则引物特异性有问题
                if sum(amplimer_counts) > 1:
                    primer_pair['specificity_warning'] = (
                        f"These primers may amplify {sum(amplimer_counts)-1} "
                        f"non-target locations in addition to the target gene."
                    )
                else:
                    primer_pair['specificity'] = "Specific to target gene"
            else:
                primer_pair['specificity'] = "No amplification found (check primer design)"

    def save_primers(self, primer_results, format='fasta'):
        """
        保存引物序列到文件
        
        参数:
            primer_results (dict): design_primers()的输出结果
            format (str): 输出格式(fasta或tsv)
        """
        gene_id = primer_results['gene_id']
        output_file = os.path.join(self.output_dir, f"{gene_id}_primers.{format}")
        
        if format == 'fasta':
            with open(output_file, 'w') as f:
                for i, pair in enumerate(primer_results['primers'], 1):
                    f.write(f">pair{i}_forward_{gene_id}\n{pair['forward']['seq']}\n")
                    f.write(f">pair{i}_reverse_{gene_id}\n{pair['reverse']['seq']}\n")
        
        elif format == 'tsv':
            with open(output_file, 'w') as f:
                f.write("Pair\tType\tSequence\tStart\tLength\tTM\tGC%\tSpecificity\n")
                for i, pair in enumerate(primer_results['primers'], 1):
                    # 正向引物
                    f.write(f"{i}\tForward\t{pair['forward']['seq']}\t")
                    f.write(f"{pair['forward']['start']}\t{pair['forward']['length']}\t")
                    f.write(f"{pair['forward']['tm']:.1f}\t{pair['forward']['gc']:.1f}\t")
                    f.write(f"{pair.get('specificity', 'Not checked')}\n")
                    
                    # 反向引物
                    f.write(f"{i}\tReverse\t{pair['reverse']['seq']}\t")
                    f.write(f"{pair['reverse']['start']}\t{pair['reverse']['length']}\t")
                    f.write(f"{pair['reverse']['tm']:.1f}\t{pair['reverse']['gc']:.1f}\t")
                    f.write(f"{pair.get('specificity', 'Not checked')}\n")
        
        return output_file
