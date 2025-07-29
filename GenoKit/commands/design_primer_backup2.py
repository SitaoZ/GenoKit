import re
import warnings
import numpy as np
from Bio.Seq import Seq
from typing import Dict
import matplotlib.pyplot as plt
from matplotlib.text import Text
from Bio.SeqUtils import GC, MeltingTemp as mt
from matplotlib.patches import Rectangle, Arrow

warnings.filterwarnings('ignore', category=UserWarning, module='Bio.SeqUtils.MeltingTemp')

class PrimerDesigner:
    def __init__(self, gene_name: str,
                 transcripts: Dict[str, str],
                 min_gc: float = 0.4,
                 max_gc: float = 0.6,
                 min_len: int = 18,
                 max_len: int = 25,
                 min_tm: int = 55,
                 max_tm: int = 65,
                 min_size: int = 100,
                 max_size: int = 300):
        """
        初始化引物设计器
        
        参数:
            gene_name: str, 基因名称
            transcripts: dict, 键为转录本名称，值为转录本序列
            min_gc: float, GC含量最小值
            max_gc: float, GC含量最大值
            min_len: int, 引物长度
            max_len: int, 引物长度
            min_tm: int, TM值
            max_tm: int, TM值
            min_size: int, 扩增片段长度
            max_size: int, 扩增片段长度
        """
        self.transcripts = transcripts
        self.gene_name = gene_name
        # 引物GC含量范围 (min, max)
        self.min_gc = min_gc
        self.max_gc = max_gc
        # 引物长度范围 (min, max)
        self.min_len = min_len
        self.max_len = max_len
        # 引物Tm值范围 (min, max)
        self.min_tm = min_tm
        self.max_tm = max_tm
        # 期望的PCR产物大小范围 (min, max)
        self.min_size = min_size
        self.max_size = max_size
          
        self.primers = {}  # 存储设计的引物
        
    def reverse_complement(self, seq):
        """生成反向互补序列"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complement[base] for base in reversed(seq)])
    
    def has_self_complement(self, primer, min_length=4):
        """检查是否存在至少min_length个连续的互补配对"""
        rc = self.reverse_complement(primer)
        n = len(primer)
        for i in range(n - min_length + 1):
            for j in range(n - min_length + 1):
                match = True
                for k in range(min_length):
                    primer_base = primer[i + k]
                    rc_base = rc[j + k]
                    if (primer_base == 'A' and rc_base != 'T') or \
                       (primer_base == 'T' and rc_base != 'A') or \
                       (primer_base == 'C' and rc_base != 'G') or \
                       (primer_base == 'G' and rc_base != 'C'):
                        match = False
                        break
                if match:
                    return True
        return False

    def validate_primer(self, primer):
        # 检查是否只包含ATCG
        if not re.match("^[ATCGatcg]+$", primer):
            # print("错误：引物包含非法碱基")
            return False
        primer = primer.upper()
        # 检查长度
        length = len(primer)
        if not (self.min_len <= length <= self.max_len):
            # print("错误：引物长度应在18-30个碱基之间")
            return False
        # 计算GC含量
        gc = primer.count('G') + primer.count('C')
        gc_percent = gc / length
        if not (self.min_gc <= gc_percent <= self.max_gc):
            # print("错误：GC含量应在40%-60%之间", gc_percent, self.min_gc, self.max_gc)
            return False
        # 计算TM值
        tm_val = mt.Tm_Wallace(primer)
        if not (self.min_tm <= tm_val <= self.max_tm):
            # print("错误：TM值应在55-65之间", tm_val, self.min_tm, self.max_tm)
            return False
        # 检查连续三个相同碱基
        for i in range(len(primer) - 2):
            if primer[i] == primer[i+1] == primer[i+2]:
                # print("错误：存在连续三个相同碱基")
                return False
        # 检查自身互补
        if self.has_self_complement(primer):
            # print("错误：引物存在自身互补序列")
            return False
        # print("引物合理")
        return True

    def design_specific_primers(self, transcript_name):
        """
        设计特定转录本的特异引物
        
        参数:
            transcript_name: str, 目标转录本名称
        返回:
            dict, 包含正向和反向引物信息
        """
        if transcript_name not in self.transcripts:
            raise ValueError(f"Transcript {transcript_name} not found in the database")
            
        target_seq = self.transcripts[transcript_name]
        seq_length = len(target_seq)

        primer_size = (self.min_len, self.max_len)
        tm = (self.min_tm, self.max_tm)
        product_size = (self.min_size, self.max_size)
        gc_content = (self.min_gc, self.max_gc)
        # 简单的引物设计策略: 在序列两端寻找合适的引物
        # 实际应用中可以使用更复杂的算法
        
        # 设计正向引物 (5'端)
        forward_primer = None
        reverse_primer = None
        for i in range(0, seq_length//2, 1):
            for f_size in range(primer_size[0], primer_size[1]+1):
                if i + f_size > seq_length:
                    continue
                f_candidate = target_seq[i : i+f_size]
                if self.validate_primer(f_candidate):
                    forward_primer = f_candidate
                    forward_start = i
                    forward_end = i + f_size
                    print("forward_end:", forward_end)
                    # 设计反向引物 (3'端)
                    for j in range(seq_length-1, seq_length//2, -1):
                        for r_size in range(primer_size[0], primer_size[1]+1):
                            if j - r_size < 0:
                                continue
                            r_candidate = target_seq[j-r_size+1 : j+1]
                            r_candidate_rc = str(Seq(r_candidate).reverse_complement())
                            if self.validate_primer(r_candidate_rc):
                                reverse_start = j - r_size + 1
                                reverse_end = j + 1
                                # 检查产物大小
                                product_len_tmp = reverse_start - forward_end
                                if not (product_size[0] <= product_len_tmp <= product_size[1]):
                                    continue
                                else:
                                    product_len = reverse_start - forward_end
                                    reverse_primer = r_candidate_rc
                                    # print(reverse_start, forward_end, product_size[0], product_len, product_size[1])
                                    if reverse_primer:
                                        break 
                        if forward_primer and reverse_primer:
                            break
                    if forward_primer and reverse_primer:
                        break
                if forward_primer and reverse_primer:
                    break
            if forward_primer and reverse_primer:
                break
        
        if not forward_primer or not reverse_primer:
            raise ValueError("Failed to design primers meeting the specified criteria")
            
        # 检查产物大小
        product_len = reverse_start - forward_end
        if not (product_size[0] <= product_len <= product_size[1]):
            print(reverse_start, forward_end, product_size[0], product_len, product_size[1])
            raise ValueError(f"Designed primers produce a product of {product_len} bp, which is outside the desired range")
            
        # 存储引物信息
        primer_info = {
            'forward': {
                'sequence': forward_primer,
                'start': forward_start,
                'end': forward_end,
                'tm': mt.Tm_Wallace(forward_primer),
                'gc': GC(forward_primer),
                'length': len(forward_primer)
            },
            'reverse': {
                'sequence': reverse_primer,
                'start': reverse_start,
                'end': reverse_end,
                'tm': mt.Tm_Wallace(reverse_primer),
                'gc': GC(reverse_primer),
                'length': len(reverse_primer)
            },
            'product_size': product_len,
            'transcript': transcript_name
        }
        
        self.primers[f"{transcript_name}_specific"] = primer_info
        return primer_info

    def visualize_primers_1(self, primer_set_name, transcript_name=None):
        """
        可视化引物在转录本上的位置
        
        参数:
            primer_set_name: str, 引物组的名称
            transcript_name: str, 可选，指定要可视化的转录本
        """
        if primer_set_name not in self.primers:
            raise ValueError(f"Primer set {primer_set_name} not found")
            
        primer_info = self.primers[primer_set_name]
        
        if transcript_name is None:
            if primer_info['transcript'] == "common":
                transcript_name = list(self.transcripts.keys())[0]
            else:
                transcript_name = primer_info['transcript']
                
        if transcript_name not in self.transcripts:
            raise ValueError(f"Transcript {transcript_name} not found")
            
        seq = self.transcripts[transcript_name]
        seq_length = len(seq)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(12, 4))
        
        # 绘制转录本
        ax.add_patch(Rectangle((0, 0.3), seq_length, 0.4, facecolor='lightblue', edgecolor='black'))
        ax.text(seq_length/2, 0.5, transcript_name, ha='center', va='center')
        
        # 绘制引物
        f_start = primer_info['forward']['start']
        f_end = primer_info['forward']['end']
        r_start = primer_info['reverse']['start']
        r_end = primer_info['reverse']['end']
        
        # 正向引物
        ax.add_patch(Rectangle((f_start, 0.7), f_end - f_start, 0.2, facecolor='green', edgecolor='black'))
        ax.text((f_start + f_end)/2, 0.8, 'F', ha='center', va='center')
        
        # 反向引物 (绘制反向互补序列)
        ax.add_patch(Rectangle((r_start, 0.1), r_end - r_start, 0.2, facecolor='red', edgecolor='black'))
        ax.text((r_start + r_end)/2, 0.2, 'R', ha='center', va='center')
        
        # 绘制PCR产物
        ax.add_patch(Rectangle((f_end, 0.3), r_start - f_end, 0.4, facecolor='yellow', edgecolor='black', alpha=0.5))
        ax.text((f_end + r_start)/2, 0.5, f"{primer_info['product_size']}bp", ha='center', va='center')
        
        # 绘制箭头表示扩增方向
        ax.arrow(f_end, 0.8, (r_start - f_end)/2, 0, head_width=0.05, head_length=10, fc='green', ec='green')
        ax.arrow(r_start, 0.2, -(r_start - f_end)/2, 0, head_width=0.05, head_length=10, fc='red', ec='red')
        
        # 设置图形属性
        ax.set_xlim(0, seq_length)
        ax.set_ylim(0, 1)
        ax.set_title(f"Primer Visualization for {transcript_name}")
        ax.set_xlabel("Position (bp)")
        ax.set_yticks([])
        ax.grid(True, axis='x', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.show()

    def visualize_primers(self, primer_set_name, transcript_name=None):
        """
        可视化引物在转录本上的位置（优化配色版）
        新增功能：
        1. 在箭头线上标注引物序列
        2. 优化整体配色方案
        """
        if primer_set_name not in self.primers:
            raise ValueError(f"Primer set {primer_set_name} not found")

        primer_info = self.primers[primer_set_name]
        
        if transcript_name is None:
            transcript_name = primer_info['transcript']
            
        if transcript_name not in self.transcripts:
            raise ValueError(f"Transcript {transcript_name} not found")
        
        seq = self.transcripts[transcript_name]
        seq_length = len(seq)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(14, 5))
        
        # 自定义配色方案
        color_palette = {
            'transcript': '#9EC8F0',
            'forward_primer': '#4DB848',  # 柔和绿色
            'reverse_primer': '#E64646',  # 珊瑚红
            'product': '#FFD700',         # 金色
            'text': '#2F4F4F',            # 深石板灰
            'arrow': '#2E8B57'            # 海绿色
        }
        
        # 绘制转录本
        ax.add_patch(Rectangle((0, 0.3), seq_length, 0.4, 
                     facecolor=color_palette['transcript'], 
                     edgecolor='#4682B4', lw=1))
        
        # 获取引物信息
        f_start = primer_info['forward']['start']
        f_end = primer_info['forward']['end']
        r_start = primer_info['reverse']['start']
        r_end = primer_info['reverse']['end']
        product_len = primer_info['product_size']
        
        # 转录本名称
        ax.text((f_end + r_start)/2, 0.55, transcript_name,
                ha='center', va='center',
                fontsize=12, color='#2F4F4F',
                fontweight='bold')
        # 绘制正向引物
        ax.add_patch(Rectangle((f_start, 0.7), f_end - f_start, 0.2,
                     facecolor=color_palette['forward_primer'], 
                     edgecolor='#2F4F4F', lw=0.8))
        ax.text((f_start + f_end)/2, 0.8, 'F', 
                ha='center', va='center', 
                fontsize=10, color='white', 
                fontweight='bold')
        
        # 绘制反向引物
        ax.add_patch(Rectangle((r_start, 0.1), r_end - r_start, 0.2,
                     facecolor=color_palette['reverse_primer'], 
                     edgecolor='#2F4F4F', lw=0.8))
        ax.text((r_start + r_end)/2, 0.2, 'R', 
                ha='center', va='center', 
                fontsize=10, color='white', 
                fontweight='bold')
        
        # 绘制PCR产物区域
        ax.add_patch(Rectangle((f_end, 0.3), r_start - f_end, 0.4,
                     facecolor=color_palette['product'], 
                     edgecolor='#DAA520', alpha=0.4, lw=1))
        ax.text((f_end + r_start)/2, 0.5, f"{product_len}bp", 
                ha='center', va='center', 
                fontsize=10, color='#2F4F4F', 
                fontstyle='italic')
        
        # 绘制箭头和引物序列标注
        arrow_params = {
            'head_width': 0.08,
            'head_length': 20,
            'fc': color_palette['arrow'],
            'ec': color_palette['arrow'],
            'lw': 1.5
        }
        
        # 正向箭头及序列标注
        f_arrow_length = (r_start - f_end) * 0.4
        f_arrow = ax.arrow(f_end, 0.8, f_arrow_length, 0, **arrow_params)
        ax.text(f_end + f_arrow_length/2, 0.85, 
                primer_info['forward']['sequence'],
                rotation=0, ha='center', va='bottom',
                fontsize=8, color=color_palette['text'],
                fontfamily='monospace', 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
        # 反向箭头及序列标注
        r_arrow_length = -(r_start - f_end) * 0.4
        r_arrow = ax.arrow(r_start, 0.2, r_arrow_length, 0, **arrow_params)
        ax.text(r_start + r_arrow_length/2, 0.15, 
                primer_info['reverse']['sequence'],
                rotation=0, ha='center', va='top',
                fontsize=8, color=color_palette['text'],
                fontfamily='monospace',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
        # 设置图形属性
        ax.set_xlim(0 - seq_length*0.02, seq_length*1.02)
        ax.set_ylim(0, 1)
        ax.set_title(f"Primer Visualization for {transcript_name}\n", 
                     fontsize=14, color='#2F4F4F')
        ax.set_xlabel("Nucleotide Position (bp)", 
                      fontsize=10, color='#2F4F4F')
        ax.set_yticks([])
        
        # 优化网格线
        ax.grid(True, axis='x', linestyle=':', 
               color='#B0C4DE', alpha=0.7)
        
        # 调整边框颜色
        for spine in ax.spines.values():
            spine.set_color('#708090')
            spine.set_linewidth(0.8)
        
        plt.tight_layout()
        plt.show()
        
    def get_primer_info(self, primer_set_name):
        """
        获取引物信息
        
        参数:
            primer_set_name: str, 引物组的名称
            
        返回:
            dict, 引物信息
        """
        if primer_set_name not in self.primers:
            raise ValueError(f"Primer set {primer_set_name} not found")
        return self.primers[primer_set_name]
        
    def print_primer_info(self, primer_set_name):
        """
        打印引物信息
        
        参数:
            primer_set_name: str, 引物组的名称
        """
        if primer_set_name not in self.primers:
            raise ValueError(f"Primer set {primer_set_name} not found")
            
        info = self.primers[primer_set_name]
        print(f"\nPrimer Set: {primer_set_name}")
        print(f"Designed for: {info['transcript']}")
        print(f"Product size: {info['product_size']} bp\n")
        
        print("Forward Primer:")
        print(f"Sequence: {info['forward']['sequence']}")
        print(f"Position: {info['forward']['start']}-{info['forward']['end']}")
        print(f"Length: {info['forward']['length']} bp")
        print(f"Tm: {info['forward']['tm']:.1f}°C")
        print(f"GC%: {info['forward']['gc']:.1f}\n")
        
        print("Reverse Primer:")
        print(f"Sequence: {info['reverse']['sequence']}")
        print(f"Position: {info['reverse']['start']}-{info['reverse']['end']}")
        print(f"Length: {info['reverse']['length']} bp")
        print(f"Tm: {info['reverse']['tm']:.1f}°C")
        print(f"GC%: {info['reverse']['gc']:.1f}\n")


if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC",
        "NM_002": "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC"
    }

    designer = PrimerDesigner("TP53", test_transcripts, min_gc = 0.4, max_gc = 0.6, min_len = 18, max_len = 25, min_tm = 55, max_tm = 65, min_size = 2000, max_size = 3000)
    designer.design_specific_primers("NM_001")
    #designer.design_specific_primers("NM_002")
    designer.visualize_primers("NM_001_specific")
    a = designer.get_primer_info("NM_001_specific")
    print(a)
    designer.print_primer_info("NM_001_specific")

