# -*- coding: utf-8 -*-
'''
filename: extract_gene.py
date: 2025/11/28 14:46
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import logging
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenoKit.database.database import gb_iter, genome_dict
from GenoKit.database.annotation_db import (
    get_gene_symbol, iter_gene_records_with_symbols, load_annotation_database,
)
from GenoKit.utils.util import parse_output, stream_parse_output
from GenoKit.utils.gtf_gff_tool import gtf_attr_to_gff, gff_attr_to_gtf


# Jbrowser2
# gene       ID=gene_id;Name=gene_name
# mRNA       ID=transcript_id;Parent=gene_id;Name=transcript_name
# exon/CDS   Parent=transcript_id

_GENE_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Chrom', 'Start', 'End', 'Strand', 'Seq', 'Length'
]
_GENE_CSV_HEADER_LINE = ','.join(_GENE_CSV_HEADER) + '\n'
_WRITE_BUFFER_SIZE = 1024 * 1024


def _write_all_genes_csv(db, genome, output_path):
    """Extract all gene sequences and write CSV in buffered batches.

    The former implementation invoked the generic output helper once for every
    gene.  For a whole human genome this means tens of thousands of tiny Python
    lists and write calls.  CSV rows contain only scalar values, so they can be
    safely assembled directly while keeping the established file layout.
    """
    buffered_rows = []
    buffered_size = 0
    with open(output_path, 'w') as out:
        out.write(_GENE_CSV_HEADER_LINE)
        for gene_id, symbol, gene in tqdm(iter_gene_records_with_symbols(db), ncols=80, total=len(db),
                                       desc='Gene processing:', colour='GREEN'):
            sequence = str(genome[gene.chr][gene.start - 1:gene.end])
            if gene.strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            row = (f'{gene_id},{symbol},{gene.chr},{gene.start},{gene.end},'
                   f'{gene.strand},{sequence},{len(sequence)}\n')
            buffered_rows.append(row)
            buffered_size += len(row)
            if buffered_size >= _WRITE_BUFFER_SIZE:
                out.write(''.join(buffered_rows))
                buffered_rows.clear()
                buffered_size = 0
        if buffered_rows:
            out.write(''.join(buffered_rows))

def get_gene(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="Gene.log"
    )
    logging.info("Loading database")
    db, t2g = load_annotation_database(args.database)
    # Fasta object
    logging.info("Loading genome")
    genome = genome_dict(args.genome)


    # print(db['ENSG00000223972'])
    if not args.gene:
        if args.output and args.output_format == 'csv':
            logging.info("Gene processing at whole-genome scale")
            _write_all_genes_csv(db, genome, args.output)
            logging.info("Gene finished at whole-genome scale")
            return

        out = open(args.output, 'w')
        if args.output_format == 'csv':
            out.write(",".join(_GENE_CSV_HEADER)+'\n')
        logging.info("Gene processing at whole-genome scale")
        for g, symbol, gene in tqdm(iter_gene_records_with_symbols(db), ncols=80, total=len(db), \
                      desc = 'Gene processing:', colour="GREEN"):
            # if isinstance(gene, dict):
            #     print(gene, g)
            gene_seq_list = []
            chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
            seq = str(genome[chrom][start-1:end])
            seq = Seq(seq)
            if gene.strand == '-':
                seq = seq.reverse_complement()

            if args.output_format == 'gtf':
                if "=" in gene.attribute:
                    # gff to gtf
                    attribute = gff_attr_to_gtf(gene.attribute)
                else:
                    attribute = gene.attribute
                gtf_format = f"{gene.chr}\t{gene.source}\t{gene.type}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.phase}\t{attribute}"
                gene_seq_list.append(gtf_format)

            elif args.output_format == 'gff':
                if "=" not in gene.attribute:
                    # gtf to gff
                    attribute = gtf_attr_to_gff(gene.attribute)
                    if 'gene_id' in attribute:
                        attribute = attribute.replace("gene_id", "ID")
                    else:
                        attribute = f"ID={g};"+attribute
                    attribute = attribute.replace("gene_name", "Name")
                else:
                    attribute = gene.attribute
                gff_format = f"{gene.chr}\t{gene.source}\t{gene.type}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.phase}\t{attribute}"
                gene_seq_list.append(gff_format)
            elif args.output_format == 'fasta':
                # fasta (defalut)
                desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                              gene.start,
                                                              gene.end,
                                                              len(seq))
                geneRecord = SeqRecord(seq, id=g, description=desc)
                gene_seq_list.append(geneRecord)
            else:
                seq = str(seq) # for polars write
                it = [g, symbol, gene.chr, gene.start, gene.end, gene.strand, seq, len(seq)]
                gene_seq_list.append(it)
                # gene_seq_list.append(dict((_GENE_CSV_HEADER[i], it[i]) for i in range(len(_GENE_CSV_HEADER))))
            stream_parse_output(args, gene_seq_list, out)
        out.close()
        logging.info("Gene finished at whole-genome scale")
    else:
        logging.info("Gene processing at gene scale")
        gene_seq_list = []
        if not db.get(args.gene):
            print('gene id not exist in database, please check.')
            sys.exit(1)
        gene = db[args.gene]['gene']
        symbol = get_gene_symbol(db, args.gene)
        chrom, start, end, strand = gene.chr, gene.start, gene.end, gene.strand
        seq = str(genome[chrom][start-1:end])
        seq = Seq(seq)
        if gene.strand == '-':
            seq = seq.reverse_complement()

        it = [args.gene, symbol, gene.chr, gene.start, gene.end, gene.strand, seq, len(seq)]
        if args.output_format == 'gtf':
            if "=" in gene.attribute:
                # gff to gtf
                attribute = gff_attr_to_gtf(gene.attribute)
            else:
                attribute = gene.attribute
            gtf_format = f"{gene.chr}\t{gene.source}\t{gene.type}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.phase}\t{attribute}"
            gene_seq_list.append(gtf_format)
        elif args.output_format == 'gff':
            if "=" not in gene.attribute:
                # gtf to gff
                attribute = gtf_attr_to_gff(gene.attribute)
                attribute = attribute.replace("gene_id", "ID")
                attribute = attribute.replace("gene_name", "Name")
            else:
                attribute = gene.attribute
            gff_format = f"{gene.chr}\t{gene.source}\t{gene.type}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.phase}\t{attribute}"
            gene_seq_list.append(gff_format)
        elif args.output_format == 'fasta':
            # fasta (defalut)
            desc='strand:%s start:%d end:%d length=%d'%(gene.strand,
                                                          gene.start,
                                                          gene.end,
                                                          len(seq))
            geneRecord = SeqRecord(seq, id=args.gene, description=desc)
            gene_seq_list.append(geneRecord)
        else:
            gene_seq_list.append(dict((_GENE_CSV_HEADER[i], it[i]) for i in range(len(_GENE_CSV_HEADER))))
        parse_output(args, gene_seq_list)
        logging.info("Gene finished at gene scale")

def get_gene_gb(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    gene = []
    for record in gb_iter(args.genbank):
        for feature in record.features:
            if feature.type == 'gene': # CDS promoter UTR 
                gene_seq = ''
                for part in feature.location.parts:
                    gene_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                gene_seq_record = SeqRecord(gene_seq, id='gene:%s'%(gene_id), 
                                            description='strand %s length %d'%(
                                            feature.strand, len(gene_seq))
                                           )
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
