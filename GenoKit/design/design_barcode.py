#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: design_barcode.py
date: 2026/6/8 09:53
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import random
import itertools
from typing import List, Tuple, Dict, Optional, Union
from collections import Counter
import numpy as np

class BarcodeDesigner:
    """DNA/RNA条形码设计器"""
    
    DNA_BASES = ['A', 'C', 'G', 'T']
    RNA_BASES = ['A', 'C', 'G', 'U']
    
    def __init__(self, barcode_type: str = 'DNA'):
        """
        初始化barcode设计器
        
        参数:
            barcode_type: 'DNA' 或 'RNA'
        """
        if barcode_type.upper() not in ['DNA', 'RNA']:
            raise ValueError("barcode_type必须是'DNA'或'RNA'")
        self.barcode_type = barcode_type.upper()
        self.bases = self.DNA_BASES if self.barcode_type == 'DNA' else self.RNA_BASES
        
    def generate_random_barcodes(
        self, 
        length: int = 8, 
        count: int = 10,
        min_gc_content: float = 0.3,
        max_gc_content: float = 0.7,
        avoid_patterns: List[str] = None,
        max_homopolymer: int = 3
    ) -> List[str]:
        """
        生成随机barcode序列
        
        参数:
            length: barcode长度
            count: 生成数量
            min_gc_content: 最小GC含量
            max_gc_content: 最大GC含量
            avoid_patterns: 需要避免的序列模式
            max_homopolymer: 最大连续相同碱基数
        
        返回:
            barcode列表
        """
        if avoid_patterns is None:
            avoid_patterns = []
            
        barcodes = []
        attempts = 0
        max_attempts = count * 100
        
        while len(barcodes) < count and attempts < max_attempts:
            attempts += 1
            # 生成随机序列
            barcode = ''.join(random.choices(self.bases, k=length))
            
            # 检查GC含量
            gc_count = barcode.count('G') + barcode.count('C')
            gc_content = gc_count / length
            if not (min_gc_content <= gc_content <= max_gc_content):
                continue
                
            # 检查连续相同碱基
            if self._has_long_homopolymer(barcode, max_homopolymer):
                continue
                
            # 检查避免的模式
            if any(pattern in barcode for pattern in avoid_patterns):
                continue
                
            # 检查序列是否已存在
            if barcode not in barcodes:
                barcodes.append(barcode)
                
        if len(barcodes) < count:
            print(f"警告: 只生成了 {len(barcodes)} 个barcode（尝试次数: {attempts}）")
            
        return barcodes
    
    def _has_long_homopolymer(self, sequence: str, max_length: int) -> bool:
        """检查是否有过长的连续相同碱基"""
        if len(sequence) == 0:
            return False
            
        current_base = sequence[0]
        current_length = 1
        
        for base in sequence[1:]:
            if base == current_base:
                current_length += 1
                if current_length > max_length:
                    return True
            else:
                current_base = base
                current_length = 1
        return False
    
    def generate_combinatorial_barcodes(
        self, 
        position_pattern: List[int] = None,
        base_options: Dict[int, List[str]] = None
    ) -> List[str]:
        """
        使用组合方法生成barcode（如Illumina测序index设计）
        
        参数:
            position_pattern: 每个位置的选项数，如[4,4,4]表示3个位置每个有4种选择
            base_options: 每个位置可选的碱基，如{0: ['A','C'], 1: ['G','T']}
        
        返回:
            所有可能的barcode组合
        """
        if position_pattern and base_options:
            raise ValueError("只能指定position_pattern或base_options之一")
            
        if position_pattern:
            # 使用位置模式生成
            positions = []
            for i, options in enumerate(position_pattern):
                positions.append(self.bases[:options])
            all_barcodes = [''.join(combo) for combo in itertools.product(*positions)]
            
        elif base_options:
            # 使用自定义的碱基选项
            positions = []
            sorted_keys = sorted(base_options.keys())
            for key in sorted_keys:
                positions.append(base_options[key])
            all_barcodes = [''.join(combo) for combo in itertools.product(*positions)]
            
        else:
            # 默认：6bp barcode，每个位置4种选择
            positions = [self.bases] * 6
            all_barcodes = [''.join(combo) for combo in itertools.product(*positions)]
            
        return all_barcodes
    
    def filter_barcodes_by_distance(
        self, 
        barcodes: List[str], 
        min_hamming_distance: int = 3,
        max_barcodes: int = None
    ) -> List[str]:
        """
        通过汉明距离筛选barcode，确保barcode之间有足够的差异
        
        参数:
            barcodes: 待筛选的barcode列表
            min_hamming_distance: 最小汉明距离
            max_barcodes: 最大返回数量
        
        返回:
            筛选后的barcode列表
        """
        if not barcodes:
            return []
            
        def hamming_distance(s1: str, s2: str) -> int:
            """计算两个序列的汉明距离"""
            if len(s1) != len(s2):
                raise ValueError("序列长度必须相同")
            return sum(c1 != c2 for c1, c2 in zip(s1, s2))
        
        # 按GC含量排序（可选，使结果更均匀）
        barcodes_sorted = sorted(barcodes, 
                                 key=lambda x: (x.count('G') + x.count('C')) / len(x))
        
        selected = [barcodes_sorted[0]]
        
        for barcode in barcodes_sorted[1:]:
            # 检查与所有已选barcode的距离
            if all(hamming_distance(barcode, s) >= min_hamming_distance for s in selected):
                selected.append(barcode)
                
            if max_barcodes and len(selected) >= max_barcodes:
                break
                
        return selected
    
    def design_dual_index_barcodes(
        self,
        i7_length: int = 8,
        i5_length: int = 8,
        pairs_count: int = 96,
        min_distance_within: int = 3,
        min_distance_between: int = 5
    ) -> Tuple[List[str], List[str]]:
        """
        设计双index测序barcode（如Illumina i7和i5）
        
        参数:
            i7_length: i7 index长度
            i5_length: i5 index长度
            pairs_count: 需要多少对
            min_distance_within: 同一个index组内的最小距离
            min_distance_between: i7和i5之间的最小距离
        
        返回:
            (i7_barcodes, i5_barcodes)
        """
        # 先生成大量候选
        i7_candidates = self.generate_random_barcodes(
            length=i7_length, 
            count=pairs_count * 10,
            min_gc_content=0.4,
            max_gc_content=0.6,
            max_homopolymer=2
        )
        
        i5_candidates = self.generate_random_barcodes(
            length=i5_length, 
            count=pairs_count * 10,
            min_gc_content=0.4,
            max_gc_content=0.6,
            max_homopolymer=2
        )
        
        # 分别筛选
        i7_selected = self.filter_barcodes_by_distance(
            i7_candidates, 
            min_hamming_distance=min_distance_within,
            max_barcodes=pairs_count
        )
        
        i5_selected = self.filter_barcodes_by_distance(
            i5_candidates, 
            min_hamming_distance=min_distance_within,
            max_barcodes=pairs_count
        )
        
        # 确保i7和i5之间也有足够差异
        # 这里简单实现：如果长度相同，检查互补性
        final_i7 = []
        final_i5 = []
        
        for i7, i5 in zip(i7_selected[:pairs_count], i5_selected[:pairs_count]):
            if i7_length == i5_length:
                # 检查是否过于相似
                matches = sum(a == b for a, b in zip(i7, i5))
                if matches <= i7_length - min_distance_between:
                    final_i7.append(i7)
                    final_i5.append(i5)
            else:
                final_i7.append(i7)
                final_i5.append(i5)
                
        return final_i7[:pairs_count], final_i5[:pairs_count]
    
    def calculate_barcode_properties(self, barcode: str) -> Dict:
        """
        计算barcode的序列属性
        
        参数:
            barcode: 序列
        
        返回:
            包含各种属性的字典
        """
        gc_count = barcode.count('G') + barcode.count('C')
        at_count = len(barcode) - gc_count
        
        # 计算GC含量
        gc_content = gc_count / len(barcode) if len(barcode) > 0 else 0
        
        # 计算熔解温度（简化的Wallace法则）
        # 每个AT对贡献2°C，每个GC对贡献4°C
        tm = 2 * at_count + 4 * gc_count
        
        # 查找最长的连续相同碱基
        max_homopolymer = 0
        current_length = 1
        for i in range(1, len(barcode)):
            if barcode[i] == barcode[i-1]:
                current_length += 1
                max_homopolymer = max(max_homopolymer, current_length)
            else:
                current_length = 1
        
        # 计算二级结构潜力（简单的回文检测）
        reverse_complement = self.get_reverse_complement(barcode)
        is_palindrome = barcode == reverse_complement
        
        return {
            'sequence': barcode,
            'length': len(barcode),
            'gc_content': gc_content,
            'gc_count': gc_count,
            'at_count': at_count,
            'tm': tm,
            'max_homopolymer': max_homopolymer,
            'is_palindrome': is_palindrome,
            'reverse_complement': reverse_complement
        }
    
    def get_reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列"""
        if self.barcode_type == 'DNA':
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        else:  # RNA
            complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def check_secondary_structure(
        self, 
        barcode: str, 
        min_stem_length: int = 4
    ) -> Dict:
        """
        检查可能的二级结构（简单实现）
        
        参数:
            barcode: 序列
            min_stem_length: 最小茎长度
        
        返回:
            二级结构信息
        """
        structure_info = {
            'has_hairpin': False,
            'stem_positions': [],
            'stem_length': 0
        }
        
        seq_len = len(barcode)
        if seq_len < min_stem_length * 2:
            return structure_info
        
        # 简化检测：寻找可能的回文区域
        for i in range(seq_len - min_stem_length * 2 + 1):
            for stem_len in range(min_stem_length, (seq_len - i) // 2 + 1):
                substr = barcode[i:i+stem_len]
                rev_comp = self.get_reverse_complement(substr)
                
                # 检查下游是否有反向互补序列
                for j in range(i + stem_len + 3, seq_len - stem_len + 1):
                    if barcode[j:j+stem_len] == rev_comp:
                        structure_info['has_hairpin'] = True
                        structure_info['stem_positions'].append((i, j, stem_len))
                        structure_info['stem_length'] = max(structure_info['stem_length'], stem_len)
        
        return structure_info


# 使用示例
def main():
    # 创建DNA barcode设计器
    dna_designer = BarcodeDesigner('DNA')
    
    print("=== 随机生成10个8bp DNA barcode ===")
    barcodes = dna_designer.generate_random_barcodes(
        length=8,
        count=10,
        min_gc_content=0.4,
        max_gc_content=0.6,
        max_homopolymer=2
    )
    
    for i, barcode in enumerate(barcodes, 1):
        props = dna_designer.calculate_barcode_properties(barcode)
        print(f"{i:2d}. {barcode} | GC: {props['gc_content']:.2f} | "
              f"Tm: {props['tm']}°C | MaxHomopolymer: {props['max_homopolymer']}")
    
    print("\n=== 通过汉明距离筛选barcode ===")
    filtered = dna_designer.filter_barcodes_by_distance(barcodes, min_hamming_distance=3)
    print(f"从 {len(barcodes)} 个筛选到 {len(filtered)} 个")
    
    print("\n=== 设计双index barcode对 ===")
    i7_barcodes, i5_barcodes = dna_designer.design_dual_index_barcodes(
        pairs_count=8,
        i7_length=8,
        i5_length=8
    )
    
    print("i7 Index\t\ti5 Index")
    print("-" * 30)
    for i7, i5 in zip(i7_barcodes, i5_barcodes):
        print(f"{i7}\t\t{i5}")
    
    print("\n=== 检查单个barcode的详细属性 ===")
    test_barcode = "ATCGATCG"
    props = dna_designer.calculate_barcode_properties(test_barcode)
    for key, value in props.items():
        print(f"{key:20}: {value}")
    
    # RNA barcode设计
    print("\n=== RNA barcode设计 ===")
    rna_designer = BarcodeDesigner('RNA')
    rna_barcodes = rna_designer.generate_random_barcodes(length=10, count=5)
    for barcode in rna_barcodes:
        print(f"RNA: {barcode}")


if __name__ == "__main__":
    main()
    # 简单使用
    # designer = BarcodeDesigner('DNA')
    # barcodes = designer.generate_random_barcodes(length=8, count=20)
    #
    # # 设计96孔板用的双index
    # i7, i5 = designer.design_dual_index_barcodes(pairs_count=96)
    #
    # # 分析单个barcode
    # print("zhusitao")
    # props = designer.calculate_barcode_properties("ATCGATCG")
    #
    # print(props)