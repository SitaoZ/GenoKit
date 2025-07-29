import re
import warnings
import numpy as np
from Bio.Seq import Seq
from typing import Dict
import matplotlib.pyplot as plt
from matplotlib.text import Text
import matplotlib.colors as mcolors
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
        设计特定转录本的特异引物（优化循环版）
        优化点：
        1. 分离正向/反向引物搜索过程
        2. 使用生成器减少内存占用
        3. 添加提前终止条件
        4. 优化参数检查逻辑
        """
        if transcript_name not in self.transcripts:
            raise ValueError(f"Transcript {transcript_name} not found")
    
        target_seq = self.transcripts[transcript_name]
        seq_length = len(target_seq)
        product_min, product_max = self.min_size, self.max_size
    
        # 生成候选正向引物
        def generate_forward_candidates():
            for start in range(0, seq_length - self.min_len + 1):
                for length in range(self.min_len, self.max_len + 1):
                    end = start + length
                    if end > seq_length:
                        continue
                    primer = target_seq[start:end]
                    if self.validate_primer(primer):
                        yield (start, end, primer)
    
        # 生成候选反向引物
        def generate_reverse_candidates():
            for end in range(seq_length, self.min_len - 1, -1):
                for length in range(self.min_len, self.max_len + 1):
                    start = end - length
                    if start < 0:
                        continue
                    primer_rc = str(Seq(target_seq[start:end]).reverse_complement())
                    if self.validate_primer(primer_rc):
                        yield (start, end, primer_rc)
    
        # 优先搜索两端区域
        search_range = min(500, seq_length//4)  # 控制搜索范围
        best_pair = None
    
        # 正向引物搜索（前1/4序列区域）
        for f_start, f_end, f_primer in generate_forward_candidates():
            if f_start > search_range:
                break  # 超出前导区域则停止
            
            # 反向引物搜索（后1/4序列区域）
            for r_start, r_end, r_primer in generate_reverse_candidates():
                if r_end < seq_length - search_range:
                    break  # 超出尾部区域则停止
    
                product_size = r_start - f_end
                if product_min <= product_size <= product_max:
                    best_pair = {
                        'forward': (f_start, f_end, f_primer),
                        'reverse': (r_start, r_end, r_primer),
                        'product_size': product_size
                    }
                    break  # 找到第一个有效对即终止
            if best_pair:
                break
    
        if not best_pair:
            # 扩大搜索范围到整个序列
            for f_start, f_end, f_primer in generate_forward_candidates():
                for r_start, r_end, r_primer in generate_reverse_candidates():
                    product_size = r_start - f_end
                    if product_min <= product_size <= product_max:
                        best_pair = {
                            'forward': (f_start, f_end, f_primer),
                            'reverse': (r_start, r_end, r_primer),
                            'product_size': product_size
                        }
                        break
                if best_pair:
                    break
    
        if not best_pair:
            raise ValueError("Failed to design primers meeting the criteria")
    
        # 构建结果结构
        forward_info = {
            'sequence': best_pair['forward'][2],
            'start': best_pair['forward'][0],
            'end': best_pair['forward'][1],
            'tm': mt.Tm_Wallace(best_pair['forward'][2]),
            'gc': GC(best_pair['forward'][2]),
            'length': len(best_pair['forward'][2])
        }
    
        reverse_info = {
            'sequence': best_pair['reverse'][2],
            'start': best_pair['reverse'][0],
            'end': best_pair['reverse'][1],
            'tm': mt.Tm_Wallace(best_pair['reverse'][2]),
            'gc': GC(best_pair['reverse'][2]),
            'length': len(best_pair['reverse'][2])
        }
    
        primer_info = {
            'forward': forward_info,
            'reverse': reverse_info,
            'product_size': best_pair['product_size'],
            'transcript': transcript_name
        }
    
        self.primers[f"{transcript_name}"] = primer_info
        return primer_info
    

    def visualize_primers_single(self, primer_set_name, transcript_name=None):
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
        
    def visualize_primers_multi(self, primer_set_names, figsize=(16, 8)):
        """
        多转录本引物可视化（统一颜色方案版）
        修改点：
        1. 所有转录本使用原始脚本的固定颜色方案
        2. 通过位置和标签区分不同转录本
        """
        import matplotlib.colors as mcolors
    
        # 使用原始脚本的颜色定义
        COLOR_PALETTE = {
            'forward': '#4DB848',  # 正向引物-绿
            'reverse': '#E64646',  # 反向引物-红
            'arrow': '#2E8B57',    # 箭头-海绿
            'product': '#FFD700'    # 产物-金
        }
    
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_ylim(0, len(primer_set_names))
        ax.set_xlim(0, 0)
    
        # 主参数配置
        transcript_height = 0.4
        primer_height = 0.15
        product_alpha = 0.3
        y_spacing = 1.0
        arrow_length_ratio = 0.4
        text_offset = 0.15
    
        for idx, set_name in enumerate(primer_set_names):
            if set_name not in self.primers:
                raise ValueError(f"Primer set {set_name} not found")
                
            primer_info = self.primers[set_name]
            transcript_name = primer_info['transcript']
            seq = self.transcripts[transcript_name]
            seq_length = len(seq)
            
            base_y = idx * y_spacing
            
            # 绘制转录本
            ax.add_patch(Rectangle((0, base_y + 0.3), seq_length, transcript_height,
                         facecolor='#9EC8F0', edgecolor='#4682B4', lw=1))
                         # facecolor='#F0F8FF', edgecolor='#4682B4', lw=1))
            ax.text(seq_length/2, base_y + 0.5, f"{transcript_name}\n({seq_length}bp)",
                    ha='center', va='center', fontsize=10, color='#2F4F4F',
                    bbox=dict(facecolor='white', alpha=0.8))  # 添加背景增强可读性
            
            # 获取引物信息
            f_start = primer_info['forward']['start']
            f_end = primer_info['forward']['end']
            r_start = primer_info['reverse']['start']
            r_end = primer_info['reverse']['end']
            product_len = primer_info['product_size']
            f_seq = primer_info['forward']['sequence']
            r_seq = primer_info['reverse']['sequence']
    
            # 绘制正向引物（统一绿色）
            ax.add_patch(Rectangle((f_start, base_y + 0.7), f_end-f_start, primer_height,
                         facecolor=COLOR_PALETTE['forward'], 
                         edgecolor='#2F4F4F', lw=0.8))
            ax.text((f_start+f_end)/2, base_y + 0.7 + primer_height/2, 'F',
                    ha='center', va='center', color='white', fontsize=9)
    
            # 绘制反向引物（统一红色）
            ax.add_patch(Rectangle((r_start, base_y + 0.1), r_end-r_start, primer_height,
                         facecolor=COLOR_PALETTE['reverse'],
                         edgecolor='#2F4F4F', lw=0.8))
            ax.text((r_start+r_end)/2, base_y + 0.1 + primer_height/2, 'R',
                    ha='center', va='center', color='white', fontsize=9)
    
            # 绘制产物区域（统一金色透明）
            ax.add_patch(Rectangle((f_end, base_y + 0.3), r_start-f_end, transcript_height,
                         facecolor=(*mcolors.to_rgb(COLOR_PALETTE['product']), product_alpha),
                         edgecolor='gray', lw=0.5))
            ax.text((f_end + r_start)/2, base_y + 0.5, f"{product_len}bp",
                    ha='center', va='center', fontsize=9, color='#2F4F4F')
    
            # 绘制箭头（统一海绿色）
            arrow_params = {
                'head_width': 0.08,
                'head_length': seq_length*0.02,
                'lw': 1.2,
                'length_includes_head': True
            }
            
            # 正向箭头及序列标注
            f_arrow_length = (r_start - f_end) * arrow_length_ratio
            ax.arrow(f_end, base_y + 0.7 + primer_height/2, 
                    f_arrow_length, 0, 
                    fc=COLOR_PALETTE['arrow'], ec=COLOR_PALETTE['arrow'], **arrow_params)
            ax.text(f_end + f_arrow_length/2, 
                    base_y + 0.7 + primer_height/2 + text_offset,
                    f_seq,
                    rotation=0, ha='center', va='bottom',
                    fontsize=8, color=COLOR_PALETTE['forward'],
                    fontfamily='monospace',
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
            # 反向箭头及序列标注
            r_arrow_length = -(r_start - f_end) * arrow_length_ratio
            ax.arrow(r_start, base_y + 0.1 + primer_height/2, 
                    r_arrow_length, 0, 
                    fc=COLOR_PALETTE['arrow'], ec=COLOR_PALETTE['arrow'], **arrow_params)
            ax.text(r_start + r_arrow_length/2,
                    base_y + 0.1 + primer_height/2 - text_offset,
                    r_seq,
                    rotation=0, ha='center', va='top',
                    fontsize=8, color=COLOR_PALETTE['reverse'],
                    fontfamily='monospace',
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
            ax.set_xlim(0, max(ax.get_xlim()[1], seq_length*1.1))
    
        # 图形装饰
        ax.set_title("Multi-Transcript Primer Visualization (Uniform Colors)\n", 
                    fontsize=14, pad=20, color='#2F4F4F')
        ax.set_xlabel("Nucleotide Position (bp)", fontsize=10)
        ax.set_yticks([i*y_spacing + 0.5 for i in range(len(primer_set_names))])
        ax.set_yticklabels([self.primers[name]['transcript'] for name in primer_set_names],
                          fontsize=10, color='#2F4F4F')
        ax.grid(True, axis='x', linestyle=':', color='#B0C4DE', alpha=0.5)
        
        # 调整边框样式
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#708090')
            spine.set_linewidth(1.2)
        
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
        
    def get_all_primer_info(self):
        return self.primers
        
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
        "NM_002": "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC",
    "NM_003": "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC",
    "NM_004": "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC",
    }

    designer = PrimerDesigner("TP53", test_transcripts, min_gc = 0.4, max_gc = 0.6, min_len = 18, max_len = 25, min_tm = 55, max_tm = 65, min_size = 700, max_size = 800)
    designer.design_specific_primers("NM_001")
    designer.design_specific_primers("NM_002")
    designer.design_specific_primers("NM_003")
    designer.design_specific_primers("NM_004")
    # designer.visualize_primers_single("NM_001_specific")
    # designer.visualize_primers_multi(["NM_001_specific", "NM_002_specific", "NM_003_specific", "NM_004_specific"])
    # designer.visualize_primers_multi(["NM_001_specific", "NM_002_specific", "NM_003_specific"])
    # a = designer.get_primer_info("NM_001_specific")
    # print(a)
    b = designer.get_all_primer_info()
    print(b)

