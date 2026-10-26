#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
filename: extract_motif.py
date: 2025/12/4 16:05
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

"""
DNA/RNA Motif Finder
支持正则表达式、简并碱基、允许错配等多种搜索方式
"""
import os
import re
import sys
import argparse
import csv
from tqdm import tqdm
from pathlib import Path
from typing import List, Dict, Tuple, Set
from collections import defaultdict


_PACKAGE_ROOT = Path(__file__).resolve().parents[1]

def _resolve_input_path(value):
    """Resolve a user path, falling back to GenoKit's packaged examples."""
    requested = Path(value).expanduser()
    if requested.is_file():
        return requested
    if not requested.is_absolute() and requested.parts and requested.parts[0] == "example":
        packaged = _PACKAGE_ROOT.joinpath(*requested.parts)
        if packaged.is_file():
            return packaged
    raise FileNotFoundError(
        "PPI input file not found: %s. Relative example/... paths are also "
        "searched under %s." % (value, _PACKAGE_ROOT / "example")
    )



class MotifFinder:
    """Motif查询器"""

    # 简并碱基字典 (IUPAC标准)
    IUPAC_CODES = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'U': 'U',
        'R': '[AG]',  # 嘌呤
        'Y': '[CTU]',  # 嘧啶 (RNA中U代替T)
        'S': '[GC]',  # 强相互作用
        'W': '[ATU]',  # 弱相互作用
        'K': '[GTU]',  # 酮基
        'M': '[AC]',  # 氨基
        'B': '[CGTU]',  # 非A
        'D': '[AGTU]',  # 非C
        'H': '[ACTU]',  # 非G
        'V': '[ACG]',  # 非T/U
        'N': '[ACGTU]',  # 任何碱基
        'X': '[ACGTU]',  # 任何碱基
        '.': '[ACGTU]',  # 任何碱基
    }

    def __init__(self, case_sensitive: bool = False, is_rna: bool = False):
        self.case_sensitive = case_sensitive
        self.is_rna = is_rna
        self.complement_map = {
            'A': 'T' if not is_rna else 'U',
            'T': 'A',
            'U': 'A',
            'C': 'G',
            'G': 'C',
            'R': 'Y',
            'Y': 'R',
            'S': 'S',
            'W': 'W',
            'K': 'M',
            'M': 'K',
            'B': 'V',
            'D': 'H',
            'H': 'D',
            'V': 'B',
            'N': 'N',
            'X': 'X',
            '.': '.'
        }

    def reverse_complement(self, sequence: str) -> str:
        """返回序列的反向互补序列"""
        seq = sequence.upper() if not self.case_sensitive else sequence
        rev_comp = ''.join(self.complement_map.get(base, base) for base in reversed(seq))
        return rev_comp

    def iupac_to_regex(self, motif: str) -> str:
        """将包含简并碱基的motif转换为正则表达式"""
        regex_parts = []
        for char in motif:
            if char in self.IUPAC_CODES:
                regex_parts.append(self.IUPAC_CODES[char])
            else:
                # 如果不是IUPAC字符，按字面量处理
                regex_parts.append(re.escape(char))

        return ''.join(regex_parts)

    def find_motifs_simple(self, sequence: str, motif: str) -> List[Tuple[int, int, str]]:
        """
        简单字符串匹配查找motif
        返回: [(起始位置, 结束位置, 匹配的序列), ...]
        """
        if not self.case_sensitive:
            sequence = sequence.upper()
            motif = motif.upper()

        results = []
        motif_len = len(motif)

        for i in range(len(sequence) - motif_len + 1):
            if sequence[i:i + motif_len] == motif:
                results.append((i, i + motif_len, sequence[i:i + motif_len]))

        return results

    def find_motifs_regex(self, sequence: str, pattern: str) -> List[Tuple[int, int, str]]:
        """
        使用正则表达式查找 motif（包含重叠匹配）
        """
        if not self.case_sensitive:
            sequence = sequence.upper()

        results = []
        # A look-ahead keeps the regex engine at every start position, so
        # overlapping alternatives such as AAAAAA are all retained.
        for match in re.finditer(f'(?=({pattern}))', sequence):
            matched_seq = match.group(1)
            start = match.start()
            results.append((start, start + len(matched_seq), matched_seq))

        return results

    def find_motifs_degenerate(self, sequence: str, motif: str) -> List[Tuple[int, int, str]]:
        """
        查找包含简并碱基的motif
        """
        regex_pattern = self.iupac_to_regex(motif)
        return self.find_motifs_regex(sequence, regex_pattern)

    def find_motifs_with_mismatches(self, sequence: str, motif: str,
                                    max_mismatches: int = 1) -> List[Tuple[int, int, str, int]]:
        """
        允许错配的motif查找
        返回: [(起始位置, 结束位置, 匹配的序列, 错配数), ...]
        """
        if not self.case_sensitive:
            sequence = sequence.upper()
            motif = motif.upper()

        results = []
        motif_len = len(motif)

        for i in range(len(sequence) - motif_len + 1):
            substring = sequence[i:i + motif_len]
            mismatches = sum(1 for a, b in zip(substring, motif) if a != b)

            if mismatches <= max_mismatches:
                results.append((i, i + motif_len, substring, mismatches))

        return results

    def find_motifs_regex_with_mismatches(
            self, sequence: str, pattern: str, max_mismatches: int):
        """Match literal regex alternatives with Hamming-distance tolerance.

        ``CTACGAAAG|AAAAAA`` is treated as two motifs.  Full regular
        expressions with quantifiers or character classes do not have a single
        fixed sequence length, so fuzzy matching for those patterns is not
        defined by this command.
        """
        alternatives = pattern.split('|')
        if not alternatives or any(
                not alternative or not re.fullmatch(r'[A-Za-z]+', alternative)
                for alternative in alternatives):
            raise ValueError(
                '--regex with --mismatches supports only literal nucleotide '
                'alternatives separated by "|", for example CTACGAAAG|AAAAAA.'
            )

        if not self.case_sensitive:
            sequence = sequence.upper()
            alternatives = [alternative.upper() for alternative in alternatives]

        # Keep one row per genomic interval, using the smallest mismatch count
        # when two alternatives describe the same interval.
        best_matches = {}
        for motif in alternatives:
            motif_length = len(motif)
            for start in range(len(sequence) - motif_length + 1):
                matched_seq = sequence[start:start + motif_length]
                mismatch_num = sum(
                    base != expected for base, expected in zip(matched_seq, motif)
                )
                if mismatch_num <= max_mismatches:
                    key = (start, start + motif_length, matched_seq)
                    previous = best_matches.get(key)
                    if previous is None or mismatch_num < previous:
                        best_matches[key] = mismatch_num

        return [
            (start, end, matched_seq, mismatch_num)
            for (start, end, matched_seq), mismatch_num in sorted(
                best_matches.items(), key=lambda item: (item[0][0], item[0][1])
            )
        ]

    def find_motifs_both_strands(self, sequence: str, motif: str,
                                 search_type: str = "both") -> Dict[str, List[Tuple[int, int, str]]]:
        """
        在两条链上查找motif
        search_type: "forward", "reverse", "both"
        """
        results = defaultdict(list)

        if search_type in ["forward", "both"]:
            # 正链
            forward_matches = self.find_motifs_degenerate(sequence, motif)
            results["forward"].extend(forward_matches)

        if search_type in ["reverse", "both"]:
            # 反链
            rev_seq = self.reverse_complement(sequence)
            reverse_matches = self.find_motifs_degenerate(rev_seq, motif)

            # 将反向链的位置转换回正链坐标
            for start, end, matched_seq in reverse_matches:
                # 转换为正向链坐标
                forward_start = len(sequence) - end
                forward_end = len(sequence) - start
                forward_seq = self.reverse_complement(matched_seq)
                results["reverse"].append((forward_start, forward_end, forward_seq))

        return results

    def find_reverse_motifs_with_mismatches(
            self, sequence: str, motif: str, max_mismatches: int):
        """Find reverse-strand matches with the same mismatch rule as forward."""
        reverse_sequence = self.reverse_complement(sequence)
        reverse_matches = self.find_motifs_with_mismatches(
            reverse_sequence, motif, max_mismatches)
        results = []
        for start, end, matched_seq, mismatch_num in reverse_matches:
            # Convert coordinates and sequence back to the forward-reference
            # coordinate system, while retaining the mismatch count.
            forward_start = len(sequence) - end
            forward_end = len(sequence) - start
            forward_seq = self.reverse_complement(matched_seq)
            results.append((forward_start, forward_end, forward_seq,
                            mismatch_num))
        return results

    def find_reverse_motifs_regex(self, sequence: str, pattern: str):
        """Find regex matches on the reverse complement and map them back."""
        reverse_sequence = self.reverse_complement(sequence)
        reverse_matches = self.find_motifs_regex(reverse_sequence, pattern)
        results = []
        for start, end, matched_seq in reverse_matches:
            forward_start = len(sequence) - end
            forward_end = len(sequence) - start
            forward_seq = self.reverse_complement(matched_seq)
            results.append((forward_start, forward_end, forward_seq))
        return results

    def find_reverse_motifs_regex_with_mismatches(
            self, sequence: str, pattern: str, max_mismatches: int):
        """Find fuzzy literal-regex matches on the reverse complement."""
        reverse_sequence = self.reverse_complement(sequence)
        reverse_matches = self.find_motifs_regex_with_mismatches(
            reverse_sequence, pattern, max_mismatches)
        results = []
        for start, end, matched_seq, mismatch_num in reverse_matches:
            forward_start = len(sequence) - end
            forward_end = len(sequence) - start
            forward_seq = self.reverse_complement(matched_seq)
            results.append((forward_start, forward_end, forward_seq,
                            mismatch_num))
        return results

    def find_all_motifs(self, sequence: str, motifs: List[str],
                        search_type: str = "both") -> Dict[str, List[Tuple[int, int, str]]]:
        """
        查找多个motif
        """
        all_results = {}

        for motif in motifs:
            motif_results = self.find_motifs_both_strands(sequence, motif, search_type)
            all_results[motif] = motif_results

        return all_results


def print_results(results: List[Tuple[int, int, str]], sequence: str,
                  motif: str, title: str = "", output: str = ""):
    """漂亮地打印搜索结果"""
    if not results:
        print(f"\n{title}: {motif} 未找到匹配", file=output)
        return

    print(f"\n{title}: 找到 {len(results)} 个匹配", file=output)
    print("-" * 80, file=output)
    print(f"Motif: {motif}", file=output)
    print(f"序列长度: {len(sequence)} bp", file=output)
    print("-" * 80, file=output)

    # print("results:", results, file=output)

    for i, (start, end, matched_seq) in enumerate(results, 1):
        # 显示上下文
        context_start = max(0, start - 10)
        context_end = min(len(sequence), end + 10)
        context = sequence[context_start:context_end]

        # 高亮显示匹配部分
        if context_start < start:
            before = context[:start - context_start]
        else:
            before = ""

        if end < context_end:
            after = context[end - context_start:]
        else:
            after = ""

        highlighted = before.lower() + matched_seq + after.lower()

        print(f"匹配 #{i}:", file=output)
        print(f"  位置: {start}-{end} (长度: {end-start} bp)", file=output)
        print(f"  序列: {matched_seq}", file=output)
        print(f"  上下文: {highlighted}", file=output)


def print_results_mismatch(results: List[Tuple[int, int, str]], sequence: str,
                  motif: str, title: str = "", output: str = ""):



    """漂亮地打印搜索结果"""
    if not results:
        print(f"{title}: 未找到匹配", file=output)
        return
    print(f"\n{title}: 找到 {len(results)} 个匹配", file=output)
    print("-" * 80, file=output)
    print(f"Motif: {motif}", file=output)
    print(f"序列长度: {len(sequence)} bp", file=output)
    print("-" * 80, file=output)

    # print("results:", results, file=output)

    for i, (start, end, matched_seq, mismatch_num) in enumerate(results, 1):
        # 显示上下文
        context_start = max(0, start - 10)
        context_end = min(len(sequence), end + 10)
        context = sequence[context_start:context_end]

        # 高亮显示匹配部分
        if context_start < start:
            before = context[:start - context_start]
        else:
            before = ""

        if end < context_end:
            after = context[end - context_start:]
        else:
            after = ""

        highlighted = before.lower() + matched_seq + after.lower()

        print(f"匹配 #{i}:", file=output)
        print(f"  位置: {start}-{end} (长度: {end-start} bp)", file=output)
        print(f"  序列: {matched_seq} (mismatch: {mismatch_num})", file=output)
        print(f"  上下文: {highlighted}", file=output)


def _match_context(sequence: str, start: int, end: int, matched_seq: str) -> str:
    """Return the context string using the same highlighting as text output."""
    context_start = max(0, start - 10)
    context_end = min(len(sequence), end + 10)
    context = sequence[context_start:context_end]
    before = context[:start - context_start] if context_start < start else ''
    after = context[end - context_start:] if end < context_end else ''
    return before.lower() + matched_seq + after.lower()


def write_results_csv(results, sequence: str, writer, fasta_id: str,
                      motif: str, strand: str = 'forward',
                      mismatches: bool = False):
    """Write motif matches as CSV rows for the ``-o`` report file."""
    for index, result in enumerate(results, 1):
        if mismatches:
            start, end, matched_seq, mismatch_num = result
            sequence_value = f'{matched_seq} (mismatch: {mismatch_num})'
        else:
            start, end, matched_seq = result
            sequence_value = matched_seq
        writer.writerow([
            fasta_id,
            motif,
            strand,
            f'match#{index}',
            f'{start}-{end}',
            f'{end - start}',
            sequence_value,
            _match_context(sequence, start, end, matched_seq),
        ])


def read_fasta_sequences(input_value: str) -> List[Tuple[str, str]]:
    """Read FASTA records while retaining each record's ID.

    A direct sequence argument remains supported and receives the stable ID
    ``input_sequence`` in CSV output.
    """
    records = []
    sequence_parts = []
    fasta_id = None
    fallback_id = os.path.basename(input_value)
    with open(input_value, 'r', encoding='utf-8') as fasta_file:
        for raw_line in fasta_file:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if fasta_id is not None:
                    records.append((fasta_id, ''.join(sequence_parts)))
                fasta_id = line[1:].split(maxsplit=1)[0] or fallback_id
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if fasta_id is not None:
        records.append((fasta_id, ''.join(sequence_parts)))
    elif sequence_parts:
        records.append((fallback_id, ''.join(sequence_parts)))
    return records


def generate_example_sequence(length: int = 100) -> str:
    """生成示例DNA序列"""
    import random
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def get_motif(args):
    '''
    Motif search
    :param args:
    args.input: input fasta file
    args.motif: sequence pattern
    args.file: motif pattern file ,one motif per line
    args.output: output file
    args.type: dna or rna
    args.case-sensitive: case sensitive
    args.both-strands: search on both chains
    args.degenerate: use IUPAC base
    args.mismatches: max mismatches allowed
    args.regex: use regex pattern
    :return:
    '''

    if not args.input:
        print("Error: input must be provided.")
        return

    if not args.motif and not args.file:
        print("Error: A motif or motif file must be provided.")
        return
    input_path = _resolve_input_path(args.input)
    fasta_records = read_fasta_sequences(input_path)

    if not fasta_records:
        print('Error: input FASTA contains no sequence.')
        return

    # 读取motif
    motifs = []
    if args.motif:
        motifs.append(args.motif)

    if args.file:
        with open(args.file, 'r') as f:
            motifs.extend(line.strip() for line in f if line.strip())

    # 创建finder
    finder = MotifFinder(case_sensitive=args.case_sensitive,
                         is_rna=(args.type == 'rna'))
    # 输出文件
    csv_writer = None
    if args.output:
        out = open(args.output, 'w', newline='', encoding='utf-8')
        csv_writer = csv.writer(out)
        # Keep the requested spelling for compatibility with existing scripts.
        csv_writer.writerow(
            ['Fasta ID', 'Motif', 'Strand', 'Match', 'Positon', 'Length', 'Seq', 'Context']
        )
    else:
        out = sys.stdout

    # 执行搜索
    if csv_writer is None:
        print("=" * 80, file=out)
        total_length = sum(len(sequence) for _, sequence in fasta_records)
        print(f"序列长度: {total_length} bp",file=out)
        motif_seqs = ",".join(motifs)
        print(f"搜索 {len(motifs)} 个motif: {motif_seqs}", file=out)
        print("=" * 80, file=out)

    for motif in tqdm(motifs, ncols = 80, total=len(motifs), \
                         desc='Motif processing:', colour="GREEN"):
        if csv_writer is None:
            print(f"\nMotif: {motif}", file=out)
        for fasta_id, sequence in fasta_records:
            if csv_writer is None and len(fasta_records) > 1:
                print(f'FASTA ID: {fasta_id}', file=out)

            if args.regex and args.mismatches > 0:
                results = finder.find_motifs_regex_with_mismatches(
                    sequence, motif, args.mismatches)
                title = f"正则表达式搜索（允许{args.mismatches}个错配）"
            elif args.mismatches > 0:
                results = finder.find_motifs_with_mismatches(sequence, motif, args.mismatches)
                title = f"允许{args.mismatches}个错配"
            elif args.regex:
                results = finder.find_motifs_regex(sequence, motif)
                title = "正则表达式搜索"
            elif args.degenerate or any(c in "RYSWKMBDHVNX." for c in motif.upper()):
                results = finder.find_motifs_degenerate(sequence, motif)
                title = "简并碱基搜索"
            else:
                results = finder.find_motifs_simple(sequence, motif)
                title = "精确匹配"

            if csv_writer is not None:
                write_results_csv(results, sequence, csv_writer, fasta_id, motif,
                                  strand='forward',
                                  mismatches=bool(args.mismatches))
            elif args.mismatches:
                print_results_mismatch(results, sequence, motif, title, out)
            else:
                print_results(results, sequence, motif, title, out)

            if args.both_strands:
                if args.regex and args.mismatches > 0:
                    both_results = {
                        'forward': results,
                        'reverse': finder.find_reverse_motifs_regex_with_mismatches(
                            sequence, motif, args.mismatches),
                    }
                elif args.mismatches > 0:
                    both_results = {
                        'forward': results,
                        'reverse': finder.find_reverse_motifs_with_mismatches(
                            sequence, motif, args.mismatches),
                    }
                elif args.regex:
                    both_results = {
                        'forward': results,
                        'reverse': finder.find_reverse_motifs_regex(
                            sequence, motif),
                    }
                else:
                    both_results = finder.find_motifs_both_strands(
                        sequence, motif, 'both')
                total_matches = 0
                if csv_writer is None:
                    print("-" * 80, file=out)
                    print(f"\nMotif: {motif}\n", file=out)
                    print(f"正负链搜索: 找到 {len([m for strand in both_results for m in both_results[strand]])} 个匹配", file=out)
                    print("-" * 80, file=out)
                for strand, matches in both_results.items():
                    if matches:
                        total_matches += len(matches)
                        strand_title = f"{strand} {title}"
                        if csv_writer is not None:
                            # The primary search above already writes forward
                            # matches.  In CSV mode, append only reverse-strand
                            # matches from the two-strand pass to avoid duplicate
                            # records when ``-b`` is used.
                            if strand == 'reverse':
                                write_results_csv(matches, sequence, csv_writer,
                                                  fasta_id, motif, strand='reverse',
                                                  mismatches=bool(args.mismatches))
                        else:
                            print_results(matches, sequence, motif, strand_title, out)
                            print("-" * 80, file=out)
                if csv_writer is None:
                    print(f"{motif} 双链总计: {total_matches} 个匹配", file=out)
    if args.output:
        out.close()

def motif_example(args):
    '''
    motif example
    :param args:
    :return:
    '''
    # 示例模式
    if args.example:

        motifs = [
            "TAGCT",  # 简单序列
            "ATG",  # 起始密码子
            "TAA|TAG|TGA",  # 终止密码子 (正则表达式)
            "GNCC",  # 简并碱基：限制性内切酶位点
            "RRAC",  # 简并碱基
        ]


        print("GenoKit motif -b -c -i zz.fasta -m ATGC -t dna ")
        return

def main():
    parser = argparse.ArgumentParser(description='DNA/RNA Motif Finder')
    parser.add_argument('-s', '--input', help='输入序列文件或直接输入序列')
    parser.add_argument('-m', '--motif', help='要查找的motif模式')
    parser.add_argument('-f', '--file', help='包含多个motif的文件（每行一个）')
    parser.add_argument('-t', '--type', choices=['dna', 'rna'], default='dna',
                        help='序列类型: dna 或 rna')
    parser.add_argument('-c', '--case-sensitive', action='store_true',
                        help='大小写敏感')
    parser.add_argument('-b', '--both-strands', action='store_true',
                        help='在两条链上搜索')
    parser.add_argument('-d', '--degenerate', action='store_true',
                        help='使用IUPAC简并碱基')
    parser.add_argument('-x', '--mismatches', type=int, default=0,
                        help='允许的最大错配数')
    parser.add_argument('-e', '--regex', action='store_true',
                        help='使用正则表达式')
    parser.add_argument('--example', action='store_true',
                        help='运行示例')

    args = parser.parse_args()

    # 示例模式
    if args.example or (not args.input and not args.motif):
        print("=" * 80)
        print("Motif Finder 示例")
        print("=" * 80)

        # 生成示例序列
        example_seq = "ATGCGATCGATCGTAGCTAGCTAGCTAGCTACGATCGATCGTAGCTAGCTAGCTAGCTAGC"
        example_seq += "TACGTAGCTAGCTAGCTAGCTAGCTACGATCGATCGATCGATCGATCGATCGATCGATC"

        motifs = [
            "TAGCT",  # 简单序列
            "ATG",  # 起始密码子
            "TAA|TAG|TGA",  # 终止密码子 (正则表达式)
            "GNCC",  # 简并碱基：限制性内切酶位点
            "RRAC",  # 简并碱基
        ]

        finder = MotifFinder(case_sensitive=False, is_rna=False)

        print(f"示例序列 ({len(example_seq)} bp):")
        print(example_seq[:80] + "..." if len(example_seq) > 80 else example_seq)
        print()

        # 测试各种搜索模式
        for motif in motifs:
            print(f"\n搜索 Motif: {motif}")
            print("-" * 60)

            if "|" in motif or args.regex:
                # 正则表达式
                results = finder.find_motifs_regex(example_seq, motif)
                print_results(results, example_seq, motif, "正则表达式搜索")
            elif any(c in "RYSWKMBDHVNX." for c in motif.upper()) or args.degenerate:
                # 简并碱基
                results = finder.find_motifs_degenerate(example_seq, motif)
                print_results(results, example_seq, motif, "简并碱基搜索")
            else:
                # 简单匹配
                results = finder.find_motifs_simple(example_seq, motif)
                print_results(results, example_seq, motif, "精确匹配")

            # 测试允许错配
            if len(motif) > 3 and "|" not in motif:
                results_mm = finder.find_motifs_with_mismatches(example_seq, motif, max_mismatches=1)
                if results_mm:
                    print(f"允许1个错配: 找到 {len(results_mm)} 个匹配")
                    for start, end, seq, mm in results_mm[:3]:  # 只显示前3个
                        print(f"  {start}-{end}: {seq} (错配: {mm})")

            # 测试双链搜索
            if args.both_strands:
                both_results = finder.find_motifs_both_strands(example_seq, motif, "both")
                for strand, matches in both_results.items():
                    if matches:
                        print(f"{strand}链: 找到 {len(matches)} 个匹配")

        return

    # 正常执行模式
    if not args.input:
        print("错误: 必须提供序列")
        parser.print_help()
        return

    if not args.motif and not args.file:
        print("错误: 必须提供motif或motif文件")
        parser.print_help()
        return

    # 读取序列
    if len(args.input) < 100 and os.path.exists(args.input):
        with open(args.input, 'r') as f:
            sequence_lines = f.readlines()
        # 移除头行和空行
        sequence = ''.join(line.strip() for line in sequence_lines
                           if not line.startswith('>') and line.strip())
    else:
        sequence = args.input

    # 读取motif
    motifs = []
    if args.motif:
        motifs.append(args.motif)
    if args.file:
        with open(args.file, 'r') as f:
            motifs.extend(line.strip() for line in f if line.strip())

    # 创建finder
    finder = MotifFinder(case_sensitive=args.case_sensitive,
                         is_rna=(args.type == 'rna'))

    # 执行搜索
    print(f"序列长度: {len(sequence)} bp")
    print(f"搜索 {len(motifs)} 个motif")
    print("=" * 80)

    for motif in motifs:
        print(f"\nMotif: {motif}")

        if args.mismatches > 0:
            results = finder.find_motifs_with_mismatches(sequence, motif, args.mismatches)
            title = f"允许{args.mismatches}个错配"
        elif args.regex:
            results = finder.find_motifs_regex(sequence, motif)
            title = "正则表达式搜索"
        elif args.degenerate or any(c in "RYSWKMBDHVNX." for c in motif.upper()):
            results = finder.find_motifs_degenerate(sequence, motif)
            title = "简并碱基搜索"
        else:
            results = finder.find_motifs_simple(sequence, motif)
            title = "精确匹配"

        print_results(results, sequence, motif, title)

        if args.both_strands:
            both_results = finder.find_motifs_both_strands(sequence, motif, "both")
            total_matches = 0
            for strand, matches in both_results.items():
                if matches:
                    print(f"{strand}链: 找到 {len(matches)} 个匹配")
                    total_matches += len(matches)
                    for start, end, seq in matches[:5]:  # 只显示前5个
                        print(f"  位置 {start}-{end}: {seq}")
                    if len(matches) > 5:
                        print(f"  ... 还有 {len(matches)-5} 个")
            print(f"双链总计: {total_matches} 个匹配")
        else:
            pass



if __name__ == "__main__":
    main()
