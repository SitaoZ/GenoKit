# -*- coding: utf-8 -*-
'''
filename: design_shrna.py
date: 2025/12/5 11:20
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

# -*- coding: utf-8 -*-
'''
filename: design_shrna.py
date: 2025/12/5 11:20
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
import logging
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from datetime import datetime
import json
import csv
from GenoKit.database.kmer_database import kmer_frequencies

# shRNA设计参数配置
SHRNA_DESIGN_PARAMS = [
    ["Parameter", "Default", "Range", "Description"],
    ["Length", "19", "19-21", "Length of targeting sequence"],
    ["Loop", "CTCGAG", "CTCGAG/TTCAAGAGA", "Loop sequence for shRNA"],
    ["GC_min", "35", "30-50", "Minimum GC content (%)"],
    ["GC_max", "55", "50-70", "Maximum GC content (%)"],
    ["Tm_min", "55", "50-65", "Minimum melting temperature (°C)"],
    ["Tm_max", "65", "55-75", "Maximum melting temperature (°C)"],
    ["Score_threshold", "60", "50-100", "Minimum overall score"],
]

# 常用loop序列
LOOP_SEQUENCES = {
    "standard": "CTCGAG",  # 常用标准loop
    "miR30": "CTCGAGAAGGTATATTGCTGTTGACAGTGAGCG",  # miR-30 backbone
    "simple": "TTCAAGAGA",  # 简单loop
}


class shRNADesigner:
    """
    shRNA设计类，用于设计针对目标基因的shRNA序列

    参数:
        gene_name (str): 基因名称
        transcripts (Dict[str, set]): 转录本字典 {transcript_id: CDS_sequence集合}
        kmer_bin (Dict[str, int]): k-mer频率字典（用于特异性评估）
        shrna_length (int): shRNA靶向序列长度(默认19)
        loop_type (str): loop序列类型(默认'standard')
        min_gc (float): 最小GC含量(默认0.35)
        max_gc (float): 最大GC含量(默认0.55)
        min_tm (float): 最小Tm值(默认55)
        max_tm (float): 最大Tm值(默认65)
    """

    def __init__(self,
                 gene_name: str,
                 transcripts: Dict[str, dict],
                 strand: str = '+',
                 kmer_bin: Dict[str, int] = None,
                 shrna_length: int = 21,
                 loop_type: str = 'standard',
                 min_gc: float = 0.35,
                 max_gc: float = 0.55,
                 min_tm: float = 55.0,
                 max_tm: float = 65.0):

        # 首先设置日志（必须在其他操作之前）
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        self.logger = logging.getLogger(__name__)

        self.gene_name = gene_name
        self.strand = strand
        self.shrna_length = shrna_length
        if loop_type in LOOP_SEQUENCES:
            self.loop_type = loop_type
            self.loop_sequence = LOOP_SEQUENCES[loop_type]
        else:
            custom_loop = loop_type.upper()
            if not custom_loop or re.search(r'[^ATCG]', custom_loop):
                raise ValueError('loop must be a named loop or an ATCG sequence.')
            self.loop_type = 'custom'
            self.loop_sequence = custom_loop
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_tm = min_tm
        self.max_tm = max_tm

        # k-mer频率字典（用于特异性评估）
        self.kmer_bin = kmer_bin if kmer_bin else {}

        # 预处理转录本（现在 logger 已经存在）
        self.transcripts = transcripts
        # self.transcripts = self._preprocess_transcripts_set(transcripts)

        # 候选池
        self.candidate_pool = defaultdict(list)

        # 黑名单模式
        self.forbidden_patterns = tuple(re.compile(pattern) for pattern in (
            r'(.)\1{3}',  # 4个连续相同碱基
            r'TTTT',  # 连续4个T
            r'AAAA',  # 连续4个A
            r'GGGG',  # 连续4个G
            r'CCCC',  # 连续4个C
        ))

        # 评分权重
        self.score_weights = {
            'gc_content': 0.15,
            'tm_stability': 0.20,
            'internal_repeat': 0.15,
            'free_energy': 0.20,
            'specificity': 0.30,
        }

    def _preprocess_transcripts_set(self, transcripts: Dict[str, set]) -> Dict[str, Dict]:
        """
        预处理转录本序列并存储正负链信息

        参数:
            transcripts: 转录本字典 {transcript_id: CDS_sequence集合}

        返回:
            处理后的转录本字典 {transcript_id: {'forward': seq, 'reverse': seq}}
        """
        processed = {}

        for tid, seqs in transcripts.items():
            # 处理集合或列表
            if isinstance(seqs, set):
                seq_list = list(seqs)
            elif isinstance(seqs, list):
                seq_list = seqs
            else:
                seq_list = [seqs]

            for seq in seq_list:
                # 清洗序列，只保留ATCG
                clean_seq = re.sub(r'[^ATCGatcg]', '', seq).upper()

                # 验证序列长度
                if len(clean_seq) < self.shrna_length + 10:
                    warning = f"转录本 {tid} 序列长度不足 {self.shrna_length + 10} bp，跳过"
                    self.logger.warning(warning)
                    continue

                # 存储正链和负链序列
                try:
                    from Bio.Seq import Seq
                    processed[tid] = {
                        'forward': clean_seq,
                        'reverse': str(Seq(clean_seq).reverse_complement())
                    }
                except ImportError:
                    # 如果没有Biopython，使用简单的反向互补
                    processed[tid] = {
                        'forward': clean_seq,
                        'reverse': self._simple_reverse_complement(clean_seq)
                    }
                break  # 每个转录本只取第一个序列

        return processed



    def _simple_reverse_complement(self, seq: str) -> str:
        """简单的反向互补序列计算（不需要Biopython）"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, 'N') for base in seq[::-1])

    def _scan_strand(self, sequence: str, strand: str) -> List[Dict]:
        """
        扫描单链上的shRNA候选位点

        参数:
            sequence: 序列
            strand: 链方向 ('+' 或 '-')

        返回:
            候选位点列表
        """
        candidates = []
        seq_len = len(sequence)

        # 滑动窗口扫描
        for i in range(seq_len - self.shrna_length + 1):
            target_seq = sequence[i:i + self.shrna_length]

            # 基础过滤
            if not self._basic_filter(target_seq):
                continue

            # 计算靶向序列的特征
            gc_content = self._calc_gc_content(target_seq)
            tm = self._calc_tm(target_seq)
            free_energy = self._calc_free_energy(target_seq)
            repeat_score = self._calc_repeat_score(target_seq)

            # 构建候选对象
            candidates.append({
                'sequence': target_seq,
                'strand': strand,
                'position': (i + 1, i + self.shrna_length),
                'gc_content': gc_content,
                'tm': tm,
                'free_energy': free_energy,
                'repeat_score': repeat_score,
            })

        return candidates

    def _basic_filter(self, seq: str) -> bool:
        """基础过滤：排除黑名单序列"""

        # 检查是否包含非ATCG字符
        if re.search(r'[^ATCG]', seq):
            return False

        # 检查黑名单模式
        for pattern in self.forbidden_patterns:
            if re.search(pattern, seq):
                return False

        # 避免以T开头（影响转录效率）
        if seq.startswith('T'):
            return False

        gc_content = self._calc_gc_content(seq)
        if not self.min_gc <= gc_content <= self.max_gc:
            return False

        tm = self._calc_tm(seq)
        if not self.min_tm <= tm <= self.max_tm:
            return False

        return True

    def _calc_gc_content(self, seq: str) -> float:
        """计算GC含量（返回0-1之间的值）"""
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq)

    def _calc_tm(self, seq: str) -> float:
        """计算熔解温度（使用Wallace规则）"""
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        tm = 4 * gc_count + 2 * at_count
        return tm

    def _calc_free_energy(self, seq: str) -> float:
        """计算自由能（基于最近邻模型）"""
        nn_params = {
            'AA': -1.0, 'AT': -0.9, 'AG': -0.8, 'AC': -0.9,
            'TA': -0.6, 'TT': -1.0, 'TG': -0.9, 'TC': -0.8,
            'GA': -0.8, 'GT': -0.9, 'GG': -1.2, 'GC': -1.5,
            'CA': -0.9, 'CT': -0.8, 'CG': -1.5, 'CC': -1.1,
        }

        total_energy = 0
        for i in range(len(seq) - 1):
            dinuc = seq[i:i + 2]
            total_energy += nn_params.get(dinuc, -0.5)

        return total_energy

    def _calc_repeat_score(self, seq: str) -> float:
        """计算内部重复评分（越高表示重复越少）"""
        penalty = 0

        # 检查2-4个碱基的重复
        for length in [2, 3, 4]:
            for i in range(len(seq) - length + 1):
                motif = seq[i:i + length]
                count = seq.count(motif)
                if count > 1:
                    penalty += (count - 1) * 5

        return max(0, 100 - penalty)

    def _get_specificity_score(self, sequence: str) -> int:
        """
        获取特异性评分（基于k-mer频率）

        参数:
            sequence: 靶向序列

        返回:
            特异性评分（频率值，越低越好）
        """
        if not self.kmer_bin:
            return 0
        return self.kmer_bin.get(sequence, 0)

    def _calculate_score(self, candidate: Dict) -> float:
        """计算综合评分（0-100）"""
        score = 0

        # GC含量评分（理想范围40-50%）
        gc = candidate['gc_content']
        if 0.40 <= gc <= 0.50:
            gc_score = 100
        elif self.min_gc <= gc <= self.max_gc:
            gc_score = 80
        else:
            gc_score = max(0, 100 - abs(gc - 0.45) * 200)
        score += gc_score * self.score_weights['gc_content']

        # Tm评分
        tm = candidate['tm']
        if self.min_tm <= tm <= self.max_tm:
            tm_score = 100
        elif tm < self.min_tm:
            tm_score = max(0, (tm / self.min_tm) * 100)
        else:
            tm_score = max(0, (self.max_tm / tm) * 100)
        score += tm_score * self.score_weights['tm_stability']

        # 内部重复评分
        repeat_score = candidate.get('repeat_score', 80)
        score += repeat_score * self.score_weights['internal_repeat']

        # 自由能评分
        energy = candidate['free_energy']
        optimal_energy = -10
        deviation = abs(energy - optimal_energy)
        energy_score = max(0, 100 - deviation * 10)
        score += energy_score * self.score_weights['free_energy']

        # 特异性评分（频率越低越好）
        spec_score = candidate.get('specificity', 0)
        if spec_score > 0:
            # 频率高的序列给予惩罚
            spec_score = max(0, 100 - spec_score * 0.1)
        else:
            spec_score = 100
        score += spec_score * self.score_weights['specificity']

        return round(min(100, max(0, score)), 2)

    def _generate_oligos(self, target_seq: str) -> Tuple[str, str]:
        """
        生成shRNA寡核苷酸序列

        参数:
            target_seq: 靶向序列

        返回:
            (sense_oligo, antisense_oligo)
        """
        # 计算反义链
        antisense_target = self._simple_reverse_complement(target_seq)

        # 构建shRNA
        sense_oligo = f"{target_seq}{self.loop_sequence}{antisense_target}"

        # 添加终止信号和酶切位点
        sense_oligo = f"{sense_oligo}TTTTTG"
        antisense_oligo = self._simple_reverse_complement(sense_oligo)
        # 添加粘性末端(AgeI:A CCGGT, EcoRI:GAATTC)
        # sense_oligo = f"AATT{sense_oligo}"
        sense_oligo = f"CCGG{sense_oligo}"
        # antisense_oligo = f"GATC{antisense_oligo}"
        antisense_oligo = f"AATT{antisense_oligo}"

        return sense_oligo, antisense_oligo

    def find_candidates(self) -> None:
        """识别所有转录本的候选shRNA"""
        if not self.transcripts:
            self.logger.warning(f"基因 {self.gene_name} 没有有效的转录本序列")
            return

        # A designer instance can be reused.  Do not accumulate the same hits
        # every time design() is called.
        self.candidate_pool.clear()
        for tid, seq in self.transcripts.items():
            for candidate in self._scan_strand(seq, self.strand):
                self._add_candidate(tid, candidate)
            # # 扫描正链
            # for candidate in self._scan_strand(seqs['forward'], '+'):
            #     self._add_candidate(tid, candidate)
            #
            # # 扫描负链
            # for candidate in self._scan_strand(seqs['reverse'], '-'):
            #     self._add_candidate(tid, candidate)

        self.logger.info(f"基因 {self.gene_name} 找到 {len(self.candidate_pool)} 个候选shRNA")

    def _add_candidate(self, tid: str, candidate: Dict) -> None:
        """存储候选shRNA并合并重复序列"""
        key = (candidate['sequence'], candidate['strand'])
        self.candidate_pool[key].append({
            'transcript': tid,
            'position': candidate['position'],
            'gc_content': candidate['gc_content'],
            'tm': candidate['tm'],
            'free_energy': candidate['free_energy'],
            'repeat_score': candidate['repeat_score'],
        })

    def filter_candidates(self,
                          require_all: bool = False,
                          require_common: bool = False,
                          *, target_mode: str | None = None) -> List[Dict]:
        """
        过滤并优化候选shRNA

        参数:
            require_all: 是否要求靶向所有转录本
            require_common: 是否要求靶向超过50%的转录本（旧 API）
            target_mode: all/any；指定时优先于旧布尔参数

        返回:
            过滤后的候选列表
        """
        if target_mode is not None:
            if target_mode not in {'all', 'any'}:
                raise ValueError("target_mode must be 'all' or 'any'")
            require_all = target_mode == 'all'
            require_common = False
        final_candidates = []

        # 确定目标转录本集合
        if not self.transcripts:
            return []
        target_transcripts = set(self.transcripts.keys())
        total_transcripts = len(target_transcripts)
        frequency_by_sequence = kmer_frequencies(
            self.kmer_bin,
            (sequence for sequence, _strand in self.candidate_pool),
        )

        for (seq, strand), hits in self.candidate_pool.items():
            # 检查靶向范围
            hit_transcripts = {h['transcript'] for h in hits}
            targets_count = len(hit_transcripts)

            if require_all and (targets_count != total_transcripts):
                continue

            if require_common and (targets_count * 2 <= total_transcripts):
                continue

            # 获取特异性评分
            specificity = frequency_by_sequence.get(seq, 0)

            # 计算平均参数
            avg_gc = sum(h['gc_content'] for h in hits) / len(hits)
            avg_tm = sum(h['tm'] for h in hits) / len(hits)
            avg_energy = sum(h['free_energy'] for h in hits) / len(hits)
            avg_repeat = sum(h['repeat_score'] for h in hits) / len(hits)

            # 构建候选字典
            candidate_dict = {
                'sequence': seq,
                'strand': strand,
                'targets': targets_count,
                'total_transcripts': total_transcripts,
                'coverage': f"{targets_count}/{total_transcripts}",
                'gc_content': round(avg_gc, 3),
                'tm': round(avg_tm, 2),
                'free_energy': round(avg_energy, 2),
                'repeat_score': round(avg_repeat, 2),
                'specificity': specificity,
                'loop': "N/A",
                'kmer_freq': specificity,
                'positions': {h['transcript']: h['position'] for h in hits},
            }

            # 计算综合评分
            candidate_dict['score'] = self._calculate_score(candidate_dict)

            # 生成oligo序列
            sense_oligo, antisense_oligo = self._generate_oligos(seq)
            candidate_dict['sense_oligo'] = sense_oligo
            candidate_dict['antisense_oligo'] = antisense_oligo
            candidate_dict['loop'] = self.loop_sequence

            final_candidates.append(candidate_dict)

        # 按评分排序
        final_candidates.sort(key=lambda x: (-x['score'], -x['targets']))

        return final_candidates

    def design(self,
               require_all: bool = False,
               require_common: bool = False,
               *, target_mode: str | None = None) -> List[Dict]:
        """
        执行完整设计流程

        参数:
            require_all: 是否要求靶向所有转录本
            require_common: 是否要求靶向超过50%的转录本（旧 API）
            target_mode: all/any；指定时优先于旧布尔参数

        返回:
            设计好的shRNA列表
        """
        if target_mode not in (None, 'all', 'any'):
            raise ValueError("target_mode must be 'all' or 'any'")
        self.find_candidates()
        return self.filter_candidates(require_all, require_common, target_mode=target_mode)

    def get_top_shrnas(self, n: int = 5, require_all: bool = False) -> List[Dict]:
        """获取前N个最佳shRNA"""
        candidates = self.design(require_all=require_all)
        if not candidates:
            return []
        return candidates[:n]

    def get_consensus_shrnas(self, n: int = 5) -> List[Tuple]:
        """获取前N个最佳通用shRNA"""
        candidates = self.design(require_all=True)
        return [(sg['sequence'], sg['strand']) for sg in candidates[:n]]

    def export_results(self,
                       candidates: List[Dict],
                       output_file: str,
                       format: str = 'csv'):
        """
        导出设计结果

        参数:
            candidates: 候选shRNA列表
            output_file: 输出文件名
            format: 输出格式 ('csv' 或 'json')
        """
        if not candidates:
            self.logger.warning("没有候选shRNA可导出")
            return

        if format == 'csv':
            self._export_to_csv(candidates, output_file)
        elif format == 'json':
            self._export_to_json(candidates, output_file)
        else:
            raise ValueError(f"不支持的格式: {format}")

    def _export_to_csv(self, candidates: List[Dict], filename: str):
        """导出为CSV格式"""
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            # 写入头部信息
            writer.writerow(['#' + '=' * 80])
            writer.writerow([f'# shRNA Design Results for {self.gene_name}'])
            writer.writerow([f'# Design Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'])
            writer.writerow([f'# Design Parameters:'])
            writer.writerow([f'#   shRNA Length: {self.shrna_length} nt'])
            writer.writerow([f'#   Loop Type: {self.loop_type} ({self.loop_sequence})'])
            writer.writerow([f'#   GC Range: {self.min_gc*100:.0f}-{self.max_gc*100:.0f}%'])
            writer.writerow([f'#   Tm Range: {self.min_tm:.0f}-{self.max_tm:.0f}°C'])
            writer.writerow([f'# Total Candidates: {len(candidates)}'])
            writer.writerow(['#' + '=' * 80])
            writer.writerow([])

            # 写入列标题
            headers = ['Rank', 'Sequence', 'Strand', 'Targets', 'GC%', 'Tm(°C)',
                       'ΔG(kcal/mol)', 'Score', 'Sense_Oligo(5→3)',
                       'Antisense_Oligo(5→3)', 'Loop']
            writer.writerow(headers)

            # 写入数据
            for idx, cand in enumerate(candidates, 1):
                writer.writerow([
                    idx,
                    cand['sequence'],
                    cand['strand'],
                    cand['targets'],
                    f"{cand['gc_content']*100:.1f}",
                    f"{cand['tm']:.1f}",
                    f"{cand['free_energy']:.2f}",
                    f"{cand['score']:.1f}",
                    cand['sense_oligo'],
                    cand['antisense_oligo'],
                    cand['loop']
                ])

        self.logger.info(f"结果已导出到: {filename}")

    def _export_to_json(self, candidates: List[Dict], filename: str):
        """导出为JSON格式"""
        export_data = {
            'gene_name': self.gene_name,
            'design_date': datetime.now().isoformat(),
            'parameters': {
                'shrna_length': self.shrna_length,
                'loop_type': self.loop_type,
                'loop_sequence': self.loop_sequence,
                'min_gc': self.min_gc,
                'max_gc': self.max_gc,
                'min_tm': self.min_tm,
                'max_tm': self.max_tm,
            },
            'total_candidates': len(candidates),
            'candidates': []
        }

        for idx, cand in enumerate(candidates, 1):
            # 转换positions中的元组为列表（JSON兼容）
            positions_json = {}
            for k, v in cand['positions'].items():
                positions_json[k] = list(v) if isinstance(v, tuple) else v

            export_data['candidates'].append({
                'rank': idx,
                'sequence': cand['sequence'],
                'strand': cand['strand'],
                'targets': cand['targets'],
                'coverage': cand['coverage'],
                'gc_content': cand['gc_content'],
                'tm': cand['tm'],
                'free_energy': cand['free_energy'],
                'score': cand['score'],
                'sense_oligo': cand['sense_oligo'],
                'antisense_oligo': cand['antisense_oligo'],
                'loop': cand['loop'],
                'positions': positions_json
            })

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2, ensure_ascii=False)

        self.logger.info(f"结果已导出到: {filename}")

    def print_summary(self, candidates: List[Dict], top_n: int = 5):
        """打印设计摘要"""
        if not candidates:
            print(f"\n基因 {self.gene_name} 没有找到合适的shRNA")
            return

        print(f"\n{'='*80}")
        print(f"shRNA设计结果 - {self.gene_name}")
        print(f"{'='*80}")
        print(f"总候选数: {len(candidates)}")
        print(f"设计参数: {self.shrna_length}nt, {self.loop_type} loop")
        print(f"GC范围: {self.min_gc*100:.0f}-{self.max_gc*100:.0f}%")
        print(f"Tm范围: {self.min_tm:.0f}-{self.max_tm:.0f}°C")
        print(f"\nTop {top_n} 推荐shRNA:")
        print(f"{'-'*80}")

        for idx, cand in enumerate(candidates[:top_n], 1):
            print(f"\n{idx}. 序列: {cand['sequence']} ({cand['strand']}链)")
            print(f"   靶向转录本: {cand['targets']}/{cand.get('total_transcripts', '?')}")
            print(f"   GC含量: {cand['gc_content']*100:.1f}%")
            print(f"   Tm: {cand['tm']:.1f}°C")
            print(f"   自由能: {cand['free_energy']:.2f} kcal/mol")
            print(f"   综合评分: {cand['score']:.1f}/100")
            print(f"   正义链: {cand['sense_oligo'][:60]}...")


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(description='shRNA设计工具')
    parser.add_argument('-i', '--input', required=True, help='输入JSON文件（包含转录本序列）')
    parser.add_argument('-g', '--gene', required=True, help='基因名称')
    parser.add_argument('-l', '--length', type=int, default=19, help='shRNA长度（默认19）')
    parser.add_argument('--loop', choices=['standard', 'miR30', 'simple'],
                        default='standard', help='Loop类型（默认standard）')
    parser.add_argument('--target-mode', choices=['all', 'any'], default='all',
                        help='all: 靶向同一基因全部输入转录本（默认）；any: 至少一个转录本')
    parser.add_argument('-n', '--number', type=int, default=10, help='输出数量（默认10）')
    parser.add_argument('-o', '--output', default='shrna_results', help='输出文件前缀')
    parser.add_argument('--format', choices=['csv', 'json'], default='csv',
                        help='输出格式（默认csv）')

    args = parser.parse_args()

    # 读取输入文件
    import json
    with open(args.input, 'r') as f:
        transcripts = json.load(f)

    # 创建设计器
    designer = shRNADesigner(
        gene_name=args.gene,
        transcripts=transcripts,
        shrna_length=args.length,
        loop_type=args.loop
    )

    # 设计shRNA
    candidates = designer.design(target_mode=args.target_mode)

    # 打印摘要
    designer.print_summary(candidates, top_n=min(args.number, 5))

    # 导出结果
    if candidates:
        output_file = f"{args.output}_{args.gene}.{args.format}"
        designer.export_results(candidates[:args.number], output_file, args.format)
        print(f"\n完整结果已保存到: {output_file}")


if __name__ == "__main__":
    # 使用示例
    test_transcripts = {
        "NM_001": {"ATGCTAGCTAGCTAGCTAGCTCGATCGATCGATCGATCGATAGCTAGCTAGC"},
        "NM_002": {"ATCGATCGATCGATCGATCGATAGCTAGCTAGCTAGCTAGC"},
        "NM_003": {"ATGAGATCCTGGCCTGGACCTGAAATGGGCACAAGGTTGTTCTTCTATGTGG"},
    }

    # 示例1：基础设计
    designer = shRNADesigner("TP53", test_transcripts, strand='+', shrna_length=19)
    candidates = designer.design()

    # 正确使用 'targets' 键名
    for shrna in candidates[:3]:
        print(f"序列: {shrna['sequence']}, 靶向转录本数: {shrna['targets']}, 评分: {shrna['score']}")

    designer.print_summary(candidates)
    designer.export_results(candidates[:5], "TP53_shrna_results.csv", "csv")

    # 示例2：要求靶向所有转录本
    universal_shrnas = designer.get_top_shrnas(n=3, require_all=True)
    print(f"\n通用shRNA: {len(universal_shrnas)}条")
    for shrna in universal_shrnas:
        print(f"  {shrna['sequence']} ({shrna['strand']}) - 靶向: {shrna['targets']}个转录本 - Score: {shrna['score']}")
