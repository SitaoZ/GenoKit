#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: extract_barcode.py
date: 2026/6/8 10:10
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import logging
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.utils.util import parse_output
from GenoKit.design.design_barcode import BarcodeDesigner


_Barcode_CSV_HEADER = [
    'Sequence', 'Length', 'GC', 'GC count', 'AT count', 'TM',
    'Max homopolymer', 'Palindrome', 'Reverse complement'
]



def get_barcode(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Barcode.log"
    )

    barcode_seq_list = []
    logging.info(f"Perform barcode design")
    nucleotide_type = args.style.upper()
    barcode_designer = BarcodeDesigner(nucleotide_type)

    barcode_count = args.count
    barcode_length = args.length
    barcodes = barcode_designer.generate_random_barcodes(
        length = barcode_length,
        count = barcode_count,
        min_gc_content = 0.4,
        max_gc_content = 0.6,
        max_homopolymer = 2
    )

    for i, barcode in tqdm(enumerate(barcodes, 1), ncols=80, total=barcode_count,
                           desc='Barcode processing:', colour="GREEN"):
        props = barcode_designer.calculate_barcode_properties(barcode)
        # print(props)
        sequence = barcode
        length = props['length']
        gc = round(props['gc_content'],2)
        gc_count = props['gc_count']
        at_count = props['at_count']
        tm = props['tm']
        max_homopolymer = props['max_homopolymer']
        is_palindrome = props['is_palindrome']
        reverse_complement = props['reverse_complement']

        if args.output_format == 'fasta':
            desc = f"length:{length}, gc:{gc}, tm:{tm}, max_homopolymer:{max_homopolymer}, is_palindrome:{is_palindrome}, reverse_complement:{reverse_complement}"
            sgRecord = SeqRecord(Seq(sequence), id='barcode'+str(i), description=desc)
            barcode_seq_list.append(sgRecord)
        else:
            it = [sequence, length, gc, gc_count, at_count, tm, max_homopolymer, is_palindrome, reverse_complement]
            barcode_seq_list.append(dict((_Barcode_CSV_HEADER[i], it[i]) for i in range(len(_Barcode_CSV_HEADER))))
    parse_output(args, barcode_seq_list)

    logging.info(f"Finished barcode design")


if __name__ == "__main__":
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
