#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: design_primer.py
date: 2025/12/12 20:35
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Arrow
from typing import Dict, List, Tuple, Optional
from typing import Optional, Tuple, Dict
import bisect
from GenoKit.database.kmer_database import kmer_frequencies


class PrimerDesigner:
    """
    极限加速优化版 PrimerDesigner
    - 取消 Biopython Tm 和 GC
    - 全局滑动窗口预计算
    - 快速发卡检测
    """

    __slots__ = (
        "transcripts", "gene_name", "min_gc", "max_gc", "min_len", "max_len",
        "min_tm", "max_tm", "min_size", "max_size", "primers",
        "kmer_bin", "kmer_size", "specificity_checker", "template_type",
        "_target_stats_cache"
    )

    def __init__(self, gene_name: str,
                 transcripts: Dict[str, str],
                 kmer_bin: Dict[str, int],
                 kmer_size: int,
                 specificity_checker=None,
                 min_gc: float = 0.4,
                 max_gc: float = 0.6,
                 min_len: int = 18,
                 max_len: int = 25,
                 min_tm: int = 55,
                 max_tm: int = 65,
                 min_size: int = 100,
                 max_size: int = 300,
                 template_type: str = "cdna"):

        if template_type not in {"genomic", "cdna"}:
            raise ValueError("template_type must be genomic or cdna")
        self.template_type = template_type
        self.transcripts = transcripts
        self.gene_name = gene_name
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_len = min_len
        self.max_len = max_len
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.min_size = min_size
        self.max_size = max_size
        self.kmer_bin = kmer_bin
        self.kmer_size = kmer_size
        self.specificity_checker = specificity_checker
        self.primers = {}
        # Only final-primer summaries are cached, never genome-wide positions.
        self._target_stats_cache = {}

    @staticmethod
    def _iter_match_positions(sequence, query):
        """Yield overlapping exact matches in an oriented template."""
        start = sequence.find(query)
        while start != -1:
            yield start
            start = sequence.find(query, start + 1)

    def _target_stats(self, primer):
        """Count both orientations, matching canonical k-mer counting.

        Each physical window is counted once (including palindromic primers).
        A template contributes at most once to coverage, even with repeats.
        The target universe is all supplied transcripts, or the genomic gene
        template; callers must keep this universe fixed during a design run.
        """
        primer = primer.upper()
        reverse = self.reverse_complement(primer)
        key = min(primer, reverse)
        if key not in self._target_stats_cache:
            occurrences = covered = 0
            for sequence in self.transcripts.values():
                sequence = sequence.upper()
                count = sum(1 for _ in self._iter_match_positions(sequence, primer))
                if reverse != primer:
                    count += sum(1 for _ in self._iter_match_positions(sequence, reverse))
                occurrences += count
                covered += bool(count)
            total = len(self.transcripts)
            self._target_stats_cache[key] = {
                'target_template_count': total,
                'target_templates_hit': covered,
                'target_coverage': covered / total if total else 0.0,
                'target_exact_matches': occurrences,
            }
        return dict(self._target_stats_cache[key])

    def _off_target_count(self, primer, total, target):
        """Subtract only when source-aware canonical metadata is available.

        Raw Jellyfish/pickle files do not expose this counting contract through
        the current adapter. Preserve their raw frequency, but report NA for
        off-target counts instead of silently assuming compatible semantics.
        Nonnegative subtraction still assumes matching reference versions.
        """
        if not self.kmer_bin:
            return 'NA', 'no_database'
        rows = getattr(self.kmer_bin, '_rows', {})
        row = rows.get(len(primer))
        if row is None:
            return 'NA', 'unknown_counting'
        expected_source = 'genome' if self.template_type == 'genomic' else 'transcriptome'
        if row['source_type'] != expected_source or not row['canonical']:
            return 'NA', 'incompatible_counting'
        if not isinstance(total, (int, np.integer)):
            return 'NA', 'missing_frequency'
        if total < target:
            return 'NA', 'count_mismatch'
        return int(total) - target, 'ok'

    def _has_target_product(self, left, right, left_length, right_length):
        """Check inward-facing nonoverlapping sites without a Cartesian product."""
        for start in left:
            lower = max(start + left_length, start + self.min_size - right_length)
            upper = start + self.max_size - right_length
            index = bisect.bisect_left(right, lower)
            if index < len(right) and right[index] <= upper:
                return True
        return False

    def annotate_pair_specificity(self, info):
        """Add basic exact-match metrics to a selected pair, without re-ranking.

        Pair coverage requires the SAME template, inward-facing sites, no
        overlap and an inclusive product-size bound. Both orientations of the
        pair are considered. This is target coverage, not off-target PCR search.
        """
        for direction in ('forward', 'reverse'):
            primer = info[direction]
            primer.update(self._target_stats(primer['sequence']))
            off_target, status = self._off_target_count(
                primer['sequence'], primer.get('kmer_freq', 'NA'),
                primer['target_exact_matches'])
            primer['off_target_exact_matches'] = off_target
            primer['specificity_status'] = status

        forward = info['forward']['sequence'].upper()
        reverse = info['reverse']['sequence'].upper()
        rc_forward = self.reverse_complement(forward)
        rc_reverse = self.reverse_complement(reverse)
        hits = 0
        for sequence in self.transcripts.values():
            sequence = sequence.upper()
            sites = {query: list(self._iter_match_positions(sequence, query))
                     for query in {forward, reverse, rc_forward, rc_reverse}}
            hits += (
                self._has_target_product(sites[forward], sites[rc_reverse],
                                         len(forward), len(reverse))
                or self._has_target_product(sites[reverse], sites[rc_forward],
                                            len(reverse), len(forward)))
        total = len(self.transcripts)
        info['target_template_count'] = total
        info['target_pair_templates_hit'] = hits
        info['target_pair_coverage'] = hits / total if total else 0.0
        return info

    # -------------------- 基础功能（加速版） --------------------

    @staticmethod
    def reverse_complement(seq: str) -> str:
        t = str.maketrans("ATCG", "TAGC")
        return seq.translate(t)[::-1]

    @staticmethod
    def calc_tm(seq: str) -> int:
        """Wallace Tm, 远快于 Biopython"""
        seq = seq.upper()
        return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))

    @staticmethod
    def calc_gc(seq: str) -> float:
        return round((seq.count("G") + seq.count("C")) / len(seq),2)

    @staticmethod
    def fast_hairpin(seq: str) -> bool:
        """
        简化 & 加速版发卡检测
        条件：长度>=10 才检查
        """
        n = len(seq)
        if n < 12:
            return False

        # 简易检查：任意4bp 与后方6bp 以上区域匹配
        for i in range(n - 10):
            sub = seq[i:i+4]
            if sub in seq[i+6:]:
                return True
        return False

    # -------------------- 引物验证 --------------------

    def validate_primer(self, primer: str) -> bool:
        """极限加速版引物验证"""

        primer = primer.upper()
        L = len(primer)

        if L < self.min_len or L > self.max_len:
            return False

        # Kmer_freq describes the complete primer.  If the catalogue does not
        # contain this exact k, do not silently substitute a shorter window.
        supports_length = getattr(self.kmer_bin, "supports_length", None)
        if self.kmer_bin and supports_length and not supports_length(L):
            return False

        if re.search(r"[^ATCG]", primer):
            return False

        # GC
        gc = self.calc_gc(primer)
        if gc < self.min_gc or gc > self.max_gc:
            return False

        # Tm
        tm = self.calc_tm(primer)
        if tm < self.min_tm or tm > self.max_tm:
            return False

        # AAA, CCC 过滤
        for i in range(L - 2):
            if primer[i] == primer[i+1] == primer[i+2]:
                return False

        # 发卡结构
        if self.fast_hairpin(primer):
            return False

        if self.specificity_checker and not self.specificity_checker.passes_primer(primer):
            return False

        return True

    def kmer_frequency(self, primer: str):
        """Return the candidate frequency without changing the primer sequence."""
        if not self.kmer_bin:
            return "NA"

        get_exact = getattr(self.kmer_bin, "get_exact", None)
        value = get_exact(primer.upper(), None) if get_exact else self.kmer_bin.get(primer.upper())
        return "NA" if value is None else int(value)

    def _filter_kmer_candidates(self, candidates):
        """Keep candidates present in the selected biological reference.

        Exact queries are performed in one batch per transcript/strand.  A
        genomic catalogue therefore removes exon-junction primers, while a
        transcriptome catalogue retains them for RT-PCR/qPCR design.
        """
        if not self.kmer_bin or not candidates:
            return candidates
        frequencies = kmer_frequencies(
            self.kmer_bin, (candidate["seq"] for candidate in candidates)
        )
        retained = []
        for candidate in candidates:
            frequency = frequencies.get(candidate["seq"], 0)
            if frequency >= 1:
                candidate["kmer_freq"] = frequency
                retained.append(candidate)
        return retained

    # -------------------- 主流程 --------------------

    def design_primers_for_transcript2(
            self, transcript_name: str,
            target_region: Optional[Tuple[int, int]] = None
            ):

        if transcript_name not in self.transcripts:
            raise ValueError(f"{transcript_name} not found")

        seq = self.transcripts[transcript_name].upper()
        N = len(seq)

        # ------------------ 目标区域处理 ------------------
        if target_region is None:
            t_start, t_end = 0, N
        else:
            t_start, t_end = target_region

        buf = min(200, N // 10)
        fwd_L = max(0, t_start - buf)
        fwd_R = min(N, t_start + buf)
        rev_L = max(0, t_end - buf)
        rev_R = min(N, t_end + buf)

        best_pair = None
        best_score = 1e9

        rc_seq = self.reverse_complement(seq)

        # ------------------ 扫描正向引物 ------------------
        for f_start in range(fwd_L, fwd_R - self.min_len + 1):
            for L in range(self.min_len, self.max_len + 1):
                f_end = f_start + L
                if f_end > fwd_R:
                    break

                f_seq = seq[f_start:f_end]
                if not self.validate_primer(f_seq):
                    continue

                # ------------------ 扫描反向引物 ------------------
                for r_end in range(rev_R, rev_L + self.min_len - 1, -1):
                    for r_L in range(self.min_len, self.max_len + 1):
                        r_start = r_end - r_L
                        if r_start < rev_L:
                            break

                        r_seq_rc = rc_seq[N - r_end:N - r_start]
                        if not self.validate_primer(r_seq_rc):
                            continue

                        product = r_start - f_end
                        if self.min_size <= product <= self.max_size:

                            f_tm = self.calc_tm(f_seq)
                            r_tm = self.calc_tm(r_seq_rc)

                            score = abs(f_tm - r_tm) * 2 + abs(self.calc_gc(f_seq) - self.calc_gc(r_seq_rc)) * 0.1

                            if score < best_score:
                                best_score = score
                                best_pair = (f_start, f_end, f_seq, r_start, r_end, r_seq_rc, product)

        if best_pair is None:
            return None

        f_start, f_end, f_seq, r_start, r_end, r_seq, product = best_pair

        info = {
            "forward": {
                "sequence": f_seq,
                "start": f_start,
                "end": f_end,
                "length": len(f_seq),
                "tm": self.calc_tm(f_seq),
                "gc": self.calc_gc(f_seq),
                "kmer_freq": self.kmer_frequency(f_seq),
            },
            "reverse": {
                "sequence": r_seq,
                "start": r_start,
                "end": r_end,
                "length": len(r_seq),
                "tm": self.calc_tm(r_seq),
                "gc": self.calc_gc(r_seq),
                "kmer_freq": self.kmer_frequency(r_seq),
            },
            "product_size": product,
            "transcript": transcript_name
        }

        self.primers[transcript_name] = info
        return info

    def design_primers_for_transcript(
            self,
            transcript_name: str,
            target_region: Optional[Tuple[int, int]] = None
    ) -> Optional[Dict]:

        if transcript_name not in self.transcripts:
            raise ValueError(
                f"{transcript_name} not found"
            )

        seq = self.transcripts[transcript_name].upper()
        N = len(seq)

        # =========================================
        # target region
        # =========================================

        if target_region is None:
            t_start, t_end = 0, N
        else:
            t_start, t_end = target_region

        buf = min(
            200,
            N // 10
        )

        fwd_L = max(
            0,
            t_start - buf
        )

        fwd_R = min(
            N,
            t_start + buf
        )

        rev_L = max(
            0,
            t_end - buf
        )

        rev_R = min(
            N,
            t_end + buf
        )

        # Without a target region, either template type may contain an
        # internal amplicon. Explicit regions retain their flanking windows.
        if target_region is None:
            fwd_L, fwd_R = 0, N
            rev_L, rev_R = 0, N

        rc_seq = self.reverse_complement(seq)

        # =====================================================
        # 1. sliding window generate forward primers
        # =====================================================

        forward_candidates = []

        for start in range(
                fwd_L,
                fwd_R - self.min_len + 1
        ):

            for L in range(
                    self.min_len,
                    self.max_len + 1
            ):

                end = start + L

                if end > fwd_R:
                    break

                primer = seq[start:end]

                if not self.validate_primer(primer):
                    continue

                forward_candidates.append(
                    {
                        "start": start,
                        "end": end,
                        "seq": primer,
                        "tm": self.calc_tm(primer),
                        "gc": self.calc_gc(primer)
                    }
                )

        forward_candidates = self._filter_kmer_candidates(forward_candidates)

        if len(forward_candidates) == 0:
            return None

        # =====================================================
        # 2. sliding window generate reverse primers
        # =====================================================

        reverse_candidates = []

        for end in range(
                rev_R,
                rev_L + self.min_len - 1,
                -1
        ):

            for L in range(
                    self.min_len,
                    self.max_len + 1
            ):

                start = end - L

                if start < rev_L:
                    break

                primer = rc_seq[
                         N - end:
                         N - start
                         ]

                if not self.validate_primer(primer):
                    continue

                reverse_candidates.append(
                    {
                        "start": start,
                        "end": end,
                        "seq": primer,
                        "tm": self.calc_tm(primer),
                        "gc": self.calc_gc(primer)
                    }
                )

        reverse_candidates = self._filter_kmer_candidates(reverse_candidates)

        if len(reverse_candidates) == 0:
            return None

        # =====================================================
        # 3. Index reverse primers by their end coordinate
        # =====================================================

        reverse_candidates.sort(
            key=lambda x: x["end"]
        )

        reverse_positions = [
            x["end"]
            for x in reverse_candidates
        ]

        # =====================================================
        # numpy arrays
        # =====================================================

        reverse_starts = np.array([x["start"] for x in reverse_candidates])
        reverse_tm = np.array(
            [
                x["tm"]
                for x in reverse_candidates
            ],
            dtype=np.float32
        )

        reverse_gc = np.array(
            [
                x["gc"]
                for x in reverse_candidates
            ],
            dtype=np.float32
        )

        # =====================================================
        # 4. binary search + numpy scoring
        # =====================================================

        best_pair = None

        best_score = np.inf

        for f in forward_candidates:

            # PCR product includes both primers: reverse.end - forward.start.
            min_r_end = max(f["start"] + self.min_size, f["end"] + self.min_len)
            max_r_end = f["start"] + self.max_size
            left = bisect.bisect_left(reverse_positions, min_r_end)
            right = bisect.bisect_right(reverse_positions, max_r_end)

            if left >= right:
                continue

            # 当前范围内所有reverse

            local_tm = reverse_tm[left:right]

            local_gc = reverse_gc[left:right]

            scores = (
                    np.abs(local_tm - f["tm"])*2
                    +
                    np.abs(local_gc - f["gc"])*0.1
            )

            # Different reverse-primer lengths can still overlap the forward primer.
            scores[reverse_starts[left:right] < f["end"]] = np.inf
            idx = np.argmin(scores)

            score = float(scores[idx])

            if score < best_score:
                r = reverse_candidates[left + idx]
                best_score = score
                best_pair = (f,r,r["end"] - f["start"])

        if best_pair is None:
            return None

        # =====================================================
        # output
        # =====================================================

        f, r, product = best_pair

        info = {

            "forward": {

                "sequence": f["seq"],
                "start": f["start"],
                "end": f["end"],
                "length": len(f["seq"]),
                "tm": f["tm"],
                "gc": f["gc"],
                "kmer_freq": f["kmer_freq"] if "kmer_freq" in f else self.kmer_frequency(f["seq"]),

            },

            "reverse": {

                "sequence": r["seq"],
                "start": r["start"],
                "end": r["end"],
                "length": len(r["seq"]),
                "tm": r["tm"],
                "gc": r["gc"],
                "kmer_freq": r["kmer_freq"] if "kmer_freq" in r else self.kmer_frequency(r["seq"]),

            },

            "product_size": product,

            "score": best_score,

            "transcript": transcript_name

        }

        self.annotate_pair_specificity(info)
        self.primers[transcript_name] = info

        return info

    def design_primers_for_all_transcripts(self) -> Dict[str, Dict]:
        """为所有转录本设计引物"""
        results = {}
        for transcript_name in self.transcripts:
            try:
                primers = self.design_primers_for_transcript(transcript_name)
                results[transcript_name] = primers
            except Exception as e:
                print(f"无法为 {transcript_name} 设计引物: {e}")
                results[transcript_name] = None
        return results

    def get_primer_info(self, transcript_name: str) -> Optional[Dict]:
        """获取指定转录本的引物信息, 字典返回"""
        primer_key = f"{transcript_name}"
        return self.primers.get(primer_key)

    def get_all_primer_info(self) -> Dict:
        """获取所有引物信息"""
        return self.primers

    def print_primer_info(self, transcript_name: str):
        """打印引物信息"""
        primer_info = self.get_primer_info(transcript_name)
        if not primer_info:
            print(f"没有找到转录本 {transcript_name} 的引物信息")
            return

        print(f"\n引物组: {transcript_name}")
        print(f"设计的转录本: {primer_info['transcript']}")
        print(f"产物大小: {primer_info['product_size']} bp\n")

        print("正向引物:")
        print(f"序列: {primer_info['forward']['sequence']}")
        print(f"位置: {primer_info['forward']['start']}-{primer_info['forward']['end']}")
        print(f"长度: {primer_info['forward']['length']} bp")
        print(f"Tm: {primer_info['forward']['tm']:.1f}°C")
        print(f"GC: {primer_info['forward']['gc']:.1f}\n")

        print("反向引物:")
        print(f"序列: {primer_info['reverse']['sequence']}")
        print(f"位置: {primer_info['reverse']['start']}-{primer_info['reverse']['end']}")
        print(f"长度: {primer_info['reverse']['length']} bp")
        print(f"Tm: {primer_info['reverse']['tm']:.1f}°C")
        print(f"GC: {primer_info['reverse']['gc']:.1f}\n")

    # 可视化方法保持不变
    def visualize_primers_single(self, primer_set_name: str, transcript_name=None):
        """可视化单个转录本的引物"""
        # 这里调用原有的可视化代码
        # ... (保持原有的可视化代码不变)
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
            'product': '#FFD700',  # 金色
            'text': '#2F4F4F',  # 深石板灰
            'arrow': '#2E8B57'  # 海绿色
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
        ax.text((f_end + r_start) / 2, 0.55, transcript_name,
                ha='center', va='center',
                fontsize=12, color='#2F4F4F',
                fontweight='bold')
        # 绘制正向引物
        ax.add_patch(Rectangle((f_start, 0.7), f_end - f_start, 0.2,
                               facecolor=color_palette['forward_primer'],
                               edgecolor='#2F4F4F', lw=0.8))
        ax.text((f_start + f_end) / 2, 0.8, 'F',
                ha='center', va='center',
                fontsize=10, color='white',
                fontweight='bold')

        # 绘制反向引物
        ax.add_patch(Rectangle((r_start, 0.1), r_end - r_start, 0.2,
                               facecolor=color_palette['reverse_primer'],
                               edgecolor='#2F4F4F', lw=0.8))
        ax.text((r_start + r_end) / 2, 0.2, 'R',
                ha='center', va='center',
                fontsize=10, color='white',
                fontweight='bold')

        # 绘制PCR产物区域
        ax.add_patch(Rectangle((f_start, 0.3), r_end - f_start, 0.4,
                               facecolor=color_palette['product'],
                               edgecolor='#DAA520', alpha=0.4, lw=1))
        ax.text((f_end + r_start) / 2, 0.5, f"{product_len}bp",
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
        ax.text(f_end + f_arrow_length / 2, 0.85,
                primer_info['forward']['sequence'],
                rotation=0, ha='center', va='bottom',
                fontsize=8, color=color_palette['text'],
                fontfamily='monospace',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        # 反向箭头及序列标注
        r_arrow_length = -(r_start - f_end) * 0.4
        r_arrow = ax.arrow(r_start, 0.2, r_arrow_length, 0, **arrow_params)
        ax.text(r_start + r_arrow_length / 2, 0.15,
                primer_info['reverse']['sequence'],
                rotation=0, ha='center', va='top',
                fontsize=8, color=color_palette['text'],
                fontfamily='monospace',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        # 设置图形属性
        ax.set_xlim(0 - seq_length * 0.02, seq_length * 1.02)
        ax.set_ylim(0, 1)
        ax.set_title(f"Primer Visualization for {transcript_name} ({self.template_type})\n",
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

    def visualize_primers_multi(self, primer_set_names: List[str]):
        """可视化多个转录本的引物"""
        # ... (保持原有的可视化代码不变)
        """
                多转录本引物可视化（统一颜色方案版）
                修改点：
                1. 所有转录本使用原始脚本的固定颜色方案
                2. 通过位置和标签区分不同转录本
                """


        # 使用原始脚本的颜色定义
        COLOR_PALETTE = {
            'forward': '#4DB848',  # 正向引物-绿
            'reverse': '#E64646',  # 反向引物-红
            'arrow': '#2E8B57',  # 箭头-海绿
            'product': '#FFD700'  # 产物-金
        }
        fig, ax = plt.subplots(figsize = (16, 8))
        ax.set_ylim(0, len(primer_set_names))
        max_sequence_length = 0

        # 主参数配置
        transcript_height = 0.4
        primer_height = 0.15
        product_alpha = 0.3
        y_spacing = 1.0
        arrow_length_ratio = 0.4
        text_offset = 0.15

        for idx, set_name in enumerate(primer_set_names):
            if set_name not in self.primers:
                logging.warning(f"Primer set '{set_name}' not found, continuing anyway.")
                # raise ValueError(f"Primer set {set_name} not found")
                # continue

            primer_info = self.primers[set_name]
            transcript_name = primer_info['transcript']
            seq = self.transcripts[transcript_name]
            seq_length = len(seq)
            max_sequence_length = max(max_sequence_length, seq_length)

            base_y = idx * y_spacing

            # 绘制转录本
            ax.add_patch(Rectangle((0, base_y + 0.3), seq_length, transcript_height,
                                   facecolor='#9EC8F0', edgecolor='#4682B4', lw=1))
            # facecolor='#F0F8FF', edgecolor='#4682B4', lw=1))
            ax.text(seq_length / 2, base_y + 0.5, f"{transcript_name} [{self.template_type}]\n({seq_length}bp)",
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

            # unipue or multiple
            if self.kmer_bin and len(f_seq) == self.kmer_size and len(r_seq) == self.kmer_size:
                # f_seq_tag = self.kmer_bin[f_seq] if self.kmer_bin.get(f_seq) else ''
                # KmerCountTable should be same as kmer size
                f_seq_tag = self.kmer_bin.get(f_seq) if self.kmer_bin.get(f_seq) else ''
                f_seq = f_seq + f"({f_seq_tag})"
                # r_seq_tag = self.kmer_bin[r_seq] if self.kmer_bin.get(r_seq) else ''
                r_seq_tag = self.kmer_bin.get(r_seq) if self.kmer_bin.get(r_seq) else ''
                r_seq = r_seq + f"({r_seq_tag})"

            # 绘制正向引物（统一绿色）
            ax.add_patch(Rectangle((f_start, base_y + 0.7), f_end - f_start, primer_height,
                                   facecolor=COLOR_PALETTE['forward'],
                                   edgecolor='#2F4F4F', lw=0.8))
            ax.text((f_start + f_end) / 2, base_y + 0.7 + primer_height / 2, 'F',
                    ha='center', va='center', color='white', fontsize=9)

            # 绘制反向引物（统一红色）
            ax.add_patch(Rectangle((r_start, base_y + 0.1), r_end - r_start, primer_height,
                                   facecolor=COLOR_PALETTE['reverse'],
                                   edgecolor='#2F4F4F', lw=0.8))
            ax.text((r_start + r_end) / 2, base_y + 0.1 + primer_height / 2, 'R',
                    ha='center', va='center', color='white', fontsize=9)

            # 绘制产物区域（统一金色透明）
            ax.add_patch(Rectangle((f_start, base_y + 0.3), r_end - f_start, transcript_height,
                                   facecolor=(*mcolors.to_rgb(COLOR_PALETTE['product']), product_alpha),
                                   edgecolor='gray', lw=0.5))
            ax.text((f_end + r_start) / 2, base_y + 0.5, f"{product_len}bp",
                    ha='center', va='center', fontsize=9, color='#2F4F4F')

            # 绘制箭头（统一海绿色）
            arrow_params = {
                'head_width': 0.08,
                'head_length': seq_length * 0.02,
                'lw': 1.2,
                'length_includes_head': True
            }

            # 正向箭头及序列标注
            f_arrow_length = (r_start - f_end) * arrow_length_ratio
            ax.arrow(f_end, base_y + 0.7 + primer_height / 2,
                     f_arrow_length, 0,
                     fc=COLOR_PALETTE['arrow'], ec=COLOR_PALETTE['arrow'], **arrow_params)
            ax.text(f_end + f_arrow_length / 2,
                    base_y + 0.7 + primer_height / 2 + text_offset,
                    f_seq,
                    rotation=0, ha='center', va='bottom',
                    fontsize=8, color=COLOR_PALETTE['forward'],
                    fontfamily='monospace',
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

            # 反向箭头及序列标注
            r_arrow_length = -(r_start - f_end) * arrow_length_ratio
            ax.arrow(r_start, base_y + 0.1 + primer_height / 2,
                     r_arrow_length, 0,
                     fc=COLOR_PALETTE['arrow'], ec=COLOR_PALETTE['arrow'], **arrow_params)
            ax.text(r_start + r_arrow_length / 2,
                    base_y + 0.1 + primer_height / 2 - text_offset,
                    r_seq,
                    rotation=0, ha='center', va='top',
                    fontsize=8, color=COLOR_PALETTE['reverse'],
                    fontfamily='monospace',
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        if max_sequence_length:
            ax.set_xlim(0, max_sequence_length * 1.1)

        # 图形装饰
        ax.set_title(f"Primer Visualization ({self.template_type})\n",
                     fontsize=14, pad=20, color='#2F4F4F')
        ax.set_xlabel("Nucleotide Position (bp)", fontsize=10)
        ax.set_yticks([i * y_spacing + 0.5 for i in range(len(primer_set_names))])
        ax.set_yticklabels([self.primers[name]['transcript'] for name in primer_set_names],
                           fontsize=10, color='#2F4F4F')
        ax.grid(True, axis='x', linestyle=':', color='#B0C4DE', alpha=0.5)

        # 调整边框样式
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#708090')
            spine.set_linewidth(1.2)

        plt.tight_layout()
        # plt.show()
        plt.savefig(f'{self.gene_name}.pdf')

# 使用示例
if __name__ == "__main__":
    # 测试序列
    test_transcripts = {
        "TP53_isoform1": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGAC",
        "TP53_isoform2": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCC",
    }

    # 创建引物设计器
    designer = PrimerDesigner(
        gene_name="TP53",
        transcripts=test_transcripts,
        kmer_bin=dict(),
        kmer_size=0,
        min_gc=0.4,
        max_gc=0.6,
        min_len=18,
        max_len=25,
        min_tm=55,
        max_tm=65,
        min_size=100,
        max_size=300
    )

    # 为所有转录本设计引物
    try:
        results = designer.design_primers_for_all_transcripts()

        # 打印结果
        print('results:', results)
        for transcript_name, primer_info in results.items():
            if primer_info:
                print(transcript_name,primer_info)
                designer.print_primer_info(transcript_name)
            else:
                print(f"\n{transcript_name}: 无法设计引物")

    except Exception as e:
        print(f"设计引物时出错: {e}")
