#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: gtf_gff_tool.py
date: 2026/8/17 15:41
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import re
from urllib.parse import quote

def gtf_attr_to_gff(attr: str) -> str:
    """
    Convert GTF attributes column to GFF3-like attributes.

    Example
    -------
    gene_id "ENSG000001"; transcript_id "ENST000001"; gene_name "ABC 1";

    ->
    gene_id=ENSG000001;transcript_id=ENST000001;gene_name=ABC%201
    """

    attr = attr.strip()
    # 支持：
    # gene_id "xxx";
    # gene_id   "xxx" ;
    # gene_id xxx;
    pattern = re.compile(
        r'([^\s;]+)\s+(?:"([^"]*)"|([^;]+))\s*;?'
    )

    result = []

    for m in pattern.finditer(attr):

        key = m.group(1)

        # quoted value
        if m.group(2) is not None:
            value = m.group(2)
        else:
            value = m.group(3).strip()

        # GFF3 中空格等特殊字符进行 URL encoding
        value = quote(
            value,
            safe="._:-|,+"
        )

        result.append(f"{key}={value}")

    return ";".join(result)

def gff_attr_to_gtf(attr: str) -> str:
    """
    Convert GFF/GFF3 attributes to GTF-style attributes.

    Example
    -------
    ID=gene1;Name=ABC1;biotype=protein_coding

    ->
    ID "gene1"; Name "ABC1"; biotype "protein_coding";
    """

    attr = attr.strip().strip(";")

    result = []

    for item in attr.split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        else:
            # 兼容非标准情况
            parts = item.split(None, 1)

            if len(parts) == 2:
                key, value = parts
            else:
                continue

        key = key.strip()
        value = value.strip()

        result.append(f'{key} "{value}";')

    return " ".join(result)