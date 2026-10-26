#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: iupac_code.py
date: 2025/12/4 12:05
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

from tabulate import tabulate


def codon():
    # codon_table.py - prints the standard genetic codon table (RNA: U instead of T)
    first_bases = ['U', 'C', 'A', 'G']
    second_bases = ['U', 'C', 'A', 'G']
    third_bases = ['U', 'C', 'A', 'G']

    aa_one_letter = [
        'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
    ]

    one_to_three = {
        'F': 'Phe', 'L': 'Leu', 'S': 'Ser', 'Y': 'Tyr', '*': 'Stop',
        'C': 'Cys', 'W': 'Trp', 'P': 'Pro', 'H': 'His', 'Q': 'Gln',
        'R': 'Arg', 'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'V': 'Val', 'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly'
    }

    codons = []
    for f in first_bases:
        for s in second_bases:
            for t in third_bases:
                codons.append(f + s + t)

    codon_to_aa3 = {codon: one_to_three[aa] for codon, aa in zip(codons, aa_one_letter)}

    def format_cell(first, second):
        items = []
        for third in third_bases:
            cod = first + second + third
            aa3 = codon_to_aa3[cod]
            items.append(f"{cod}={aa3}")
        return "\n".join(items)

    col_width = 12  # 24 # 36
    sep = " | "
    header_cells = [f"2nd={s}".center(col_width) for s in second_bases]
    print(" " * 8 + sep.join(header_cells))
    print("-" * (8 + (col_width + len(sep)) * len(second_bases) - len(sep)))

    for first in first_bases:
        cell_texts = [format_cell(first, s).splitlines() for s in second_bases]
        max_lines = max(len(c) for c in cell_texts)
        print(f"1st={first}".ljust(8), end=sep)
        for i in range(max_lines):
            row_parts = []
            for cell in cell_texts:
                part = cell[i] if i < len(cell) else ""
                row_parts.append(part.ljust(col_width))
            print(sep.join(row_parts))
            if i < max_lines - 1:
                print(" " * 8, end=sep)
        print("-" * (8 + (col_width + len(sep)) * len(second_bases) - len(sep)))


def get_iupac(args):
    '''
    IUPAC code
    :param args:
    :return:
    '''
    # https://www.bioinformatics.org/sms/iupac.html
    nucleotide = [
        ["IUPAC nucleotide code", "Base"],
        ["A", "Adenine"],
        ["C", "Cytosine"],
        ["G", "Guanine"],
        ["T", "Thymine"],
        ["U", "Uracil"],
        ["R", "A or G"],
        ["Y", "C or T"],
        ["S", "G or C"],
        ["W", "A or T"],
        ["K", "G or T"],
        ["M", "A or C"],
        ["B", "C or G or T"],
        ["D", "A or G or T"],
        ["H", "A or C or T"],
        ["V", "A or C or G"],
        ["N", "any base"],
        [". or -", "gap"]
    ]
    # https://worldwide.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/
    # https://zh.wikipedia.org/wiki/%E6%B0%A8%E5%9F%BA%E9%85%B8
    amino = [
        ["IUPAC amino acid code", "Three letter code", "Amino acid", "Molecular Weight", "Chinese name"],
        ["A", "Ala", "Alanine", "89.09", "丙氨酸"],
        ["C", "Cys", "Cysteine", "121.16","半胱氨酸"],
        ["D", "Asp", "Aspartic Acid", "133.10", "天冬氨酸"],
        ["E", "Glu", "Glutamic Acid", "147.13", "谷氨酸"],
        ["F", "Phe", "Phenylalanine", "165.19", "苯丙氨酸"],
        ["G", "Gly", "Glycine", "75.07", "甘氨酸"],
        ["H", "His", "Histidine", "155.16", "组氨酸"],
        ["I", "Ile", "Isoleucine", "131.17", "异亮氨酸"],
        ["K", "Lys", "Lysine", "146.19", "赖氨酸"],
        ["L", "Leu", "Leucine", "131.17", "亮氨酸"],
        ["M", "Met", "Methionine", "149.21", "甲硫氨酸"],
        ["N", "Asn", "Asparagine", "132.12", "天冬酰胺"],
        ["P", "Pro", "Proline", "115.13", "脯氨酸"],
        ["Q", "Gln", "Glutamine", "146.15", "谷氨酰胺"],
        ["R", "Arg", "Arginine", "174.20", "精氨酸"],
        ["S", "Ser", "Serine", "105.09", "丝氨酸"],
        ["T", "Thr", "Threonine", "119.12", "苏氨酸"],
        ["V", "Val", "Valine", "117.15", "缬氨酸"],
        ["W", "Trp", "Tryptophan", "204.23", "色氨酸"],
        ["Y", "Tyr", "Tyrosine", "181.19", "酪氨酸"]
        ]


    if args.nucleotide:
        print(tabulate(nucleotide, headers="firstrow", tablefmt="grid"))
    elif args.amino:
        print(tabulate(
            amino,
            headers="firstrow",
            tablefmt="grid",
            colalign=("left", "left", "left", "left", "right"),
        ))
    elif args.codon:
        codon()
    else:
        args.command_parser.print_help()


