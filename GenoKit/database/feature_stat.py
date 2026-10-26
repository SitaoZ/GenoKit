# -*- coding: utf-8 -*-
'''
filename: feature_stat.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import os
from GenoKit.database.database import genome_dict
from GenoKit.utils.util import utr3_type, utr5_type
from GenoKit.database.annotation_db import load_annotation_database



def _fasta_statistics(genome_path):
    """Return sequence count and total bases using the FASTA index when available."""
    fai_path = genome_path + '.fai'
    if not os.path.exists(fai_path):
        # Keep the existing behaviour for an unindexed FASTA: create its index
        # through the project's standard helper, then read the small index file.
        genome_dict(genome_path)

    chromosome_count = 0
    genome_size = 0
    with open(fai_path, 'r') as fai_file:
        for line in fai_file:
            if not line.strip():
                continue
            fields = line.split('\t', 2)
            if len(fields) < 2:
                raise ValueError(f'Invalid FASTA index: {fai_path}')
            chromosome_count += 1
            genome_size += int(fields[1])
    return chromosome_count, genome_size


def calculate_statistics(db, style):
    """Calculate all annotation counts in one pass over the database."""
    if hasattr(db, 'statistics'):
        # SQLite can calculate these values with indexed SQL aggregates and
        # does not need to materialize every transcript and feature record.
        return db.statistics(style)
    utr3_t = utr3_type(style)
    utr5_t = utr5_type(style)
    transcript_nums = 0
    coding_gene_nums = 0
    coding_transcript_nums = 0
    utr5_transcript_nums = 0
    utr3_transcript_nums = 0
    multiple_exon_transcript_nums = 0

    for gene in db.values():
        gene_has_cds = False
        for transcript_id, transcript in gene.items():
            if transcript_id in {'gene', 'symbol'}:
                continue
            transcript_nums += 1
            if transcript.get('CDS'):
                coding_transcript_nums += 1
                gene_has_cds = True
            if transcript.get(utr5_t):
                utr5_transcript_nums += 1
            if transcript.get(utr3_t):
                utr3_transcript_nums += 1
            if len(transcript.get('exon', ())) >= 2:
                multiple_exon_transcript_nums += 1
        if gene_has_cds:
            coding_gene_nums += 1

    return {
        'gene_nums': len(db),
        'transcript_nums': transcript_nums,
        'coding_gene_nums': coding_gene_nums,
        'coding_transcript_nums': coding_transcript_nums,
        'utr5_transcript_nums': utr5_transcript_nums,
        'utr3_transcript_nums': utr3_transcript_nums,
        'multiple_exon_transcript_nums': multiple_exon_transcript_nums,
    }


def get_stat(args):
    '''
    parameters:
        args: parse from argparse
        args.database database created from create command
        args.genome genome fasta file
        args.output output path, if not stdout
        args.style gff of gtf
    return:
        a genome database statistics
    '''
    db, t2g = load_annotation_database(args.database)
    chrom_size, genome_size = _fasta_statistics(args.genome)
    stat = calculate_statistics(db, args.style)
    gene_nums = stat['gene_nums']
    transcript_nums = stat['transcript_nums']
    coding_gene_nums = stat['coding_gene_nums']
    coding_transcript_nums = stat['coding_transcript_nums']
    utr5_transcript_nums = stat['utr5_transcript_nums']
    utr3_transcript_nums = stat['utr3_transcript_nums']
    multiple_exon_transcript_nums = stat['multiple_exon_transcript_nums']
    if not args.output:
        print('Genome      size:', genome_size)
        print('Chromosome  size:', chrom_size)
        print('Gene       count:', gene_nums)
        print('Transcript count:', transcript_nums)
        print('Gene with CDS:', coding_gene_nums)
        print('Transcript with CDS:', coding_transcript_nums)
        print('Transcript with UTR5:', utr5_transcript_nums)
        print('Transcript with UTR3:', utr3_transcript_nums)
        print('Transcript with multiple exon:', multiple_exon_transcript_nums)
    else:
        with open(args.output, 'w') as f:
            f.write(f'Type,Count\n')
            f.write(f'Genome size,{genome_size}\n')
            f.write(f'Chromosome size,{chrom_size}\n')
            f.write(f'Gene count,{gene_nums}\n')
            f.write(f'Transcript count,{transcript_nums}\n')
            f.write(f'Gene with CDS,{coding_gene_nums}\n')
            f.write(f'Transcript with CDS,{coding_transcript_nums}\n')
            f.write(f'Transcript with UTR5,{utr5_transcript_nums}\n')
            f.write(f'Transcript with UTR3,{utr3_transcript_nums}\n')
            f.write(f'Transcript with multiple exon,{multiple_exon_transcript_nums}\n')
