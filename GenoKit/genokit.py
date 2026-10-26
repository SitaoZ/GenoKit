# -*- coding: utf-8 -*-
'''
filename: genokit.py
date: 2024/11/20 15:40
author: Sitao Zhu
mail: zhusitao1990@163.com
'''

import sys
import argparse
import importlib



# Command implementations have scientific and plotting dependencies.  Importing
# all of them just to print help makes the command-line interface needlessly slow.
_COMMAND_HANDLERS = {
    'create': ('GenoKit.database.build_db', 'create_db'),
    'fasta': ('GenoKit.database.fasta_split', 'get_fasta'),
    'iupac': ('GenoKit.database.iupac_code', 'get_iupac'),
    'kmer': ('GenoKit.database.kmer_stat', 'get_kmer'),
    'stat': ('GenoKit.database.feature_stat', 'get_stat'),

    'utr': ('GenoKit.extract.extract_utr', 'get_utr'),
    'uorf': ('GenoKit.extract.extract_uorf', 'get_uorf'),
    'cds': ('GenoKit.extract.extract_cds', 'get_cds'),
    'dorf': ('GenoKit.extract.extract_dorf', 'get_dorf'),
    'exon': ('GenoKit.extract.extract_exon', 'get_exon'),
    'intron': ('GenoKit.extract.extract_intron', 'get_intron'),
    'promoter': ('GenoKit.extract.extract_promoter', 'get_promoter'),
    'terminator': ('GenoKit.extract.extract_terminator', 'get_terminator'),
    'gene': ('GenoKit.extract.extract_gene', 'get_gene'),
    'mrna': ('GenoKit.extract.extract_mrna', 'get_mrna'),
    'transcript': ('GenoKit.extract.extract_transcript', 'get_transcript'),
    'igr': ('GenoKit.extract.extract_igr', 'get_igr'),

    'motif': ('GenoKit.design.extract_motif', 'get_motif'),
    'sirna': ('GenoKit.design.extract_sirna', 'get_sirna'),
    'sgrna': ('GenoKit.design.extract_sgrna', 'get_sgrna'),
    'shrna': ('GenoKit.design.extract_shrna', 'get_shrna'),
    'barcode': ('GenoKit.design.extract_barcode', 'get_barcode'),
    'primer': ('GenoKit.design.extract_primer', 'get_primer'),

    'circos': ('GenoKit.visualize.extract_circos', 'get_circos'),
    'hgvs': ('GenoKit.visualize.extract_hgvs', 'get_hgvs'),
    'ppi': ('GenoKit.visualize.extract_ppi', 'get_ppi'),
    'view': ('GenoKit.visualize.extract_view', 'get_view'),
}


def _lazy_handler(command):
    def handler(args):
        module_name, function_name = _COMMAND_HANDLERS[command]
        getattr(importlib.import_module(module_name), function_name)(args)

    return handler


for _command in _COMMAND_HANDLERS:
    globals()[_command] = _lazy_handler(_command)

parser = argparse.ArgumentParser()
# 添加全局参数
VERSION = "1.0.0"
parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'GenoKit v{VERSION}',
        help='Show version information'
    )


# 添加子命令
subparsers = parser.add_subparsers(help='sub-command help')


# create subcommand
create_example_text = '''
example:
  # sqlite (human.gtf.splite)
  GenoKit create -s gtf -a human.gtf -g human.fa -o ./ -p human
  GenoKit create -s gtf -a human.gtf.gz -g human.fa.gz -o ./ -p human
  GenoKit create -s gtf -a human.gtf -g human.fa -o ./ -p human --database-format sqlite
  GenoKit create -s gff -a yeast.gff -g yeast.fa -o ./ -p yeast --database-format sqlite
  # pickle (human.gtf.pkl)
  GenoKit create -s gtf -a human.gtf -g human.fa -o ./ -p human --database-format pickle
  
note: 
  • annotation file (-a) accepts either standard gzip or uncompressed format
  • reference genome (-g) requires bgzip compression for indexing support or uncompressed format
  • SQLite is the default database; legacy pickle output remains available for compatibility'''

parser_create = subparsers.add_parser('create', help='Build a genome annotation database',
                                      epilog=create_example_text,
                                      formatter_class=argparse.RawTextHelpFormatter
                                      )
parser_create.add_argument('-a', '--annotation', type=str, required=True,
                           help='genome annotation path (gff or gtf), uncompressed or gzip (.gz)')
parser_create.add_argument('-g', '--genome', type=str, required=True,
                           help='genome fasta path, uncompressed or bgzip (.gz) ONLY (for efficient indexing)')

parser_create.add_argument('-o', '--output_dir', type=str, required=True,
                           help='database output directory path')
parser_create.add_argument('-p', '--prefix', type=str, required=True,
                           help='database prefix')
parser_create.add_argument('-s', '--style', choices=['gff','gtf'],
                           required=True,
                           help='annotation syntax; use gff for GFF3 and gtf for GTF')
parser_create.add_argument('--database-format', choices=['sqlite', 'pickle'],
                           default='sqlite',
                           help='database backend (default: sqlite; pickle writes a .pkl file)')
parser_create.add_argument('--skip-fasta-index', action='store_true',
                           help='do not create a missing FASTA .fai index')
parser_create.set_defaults(func=create)


# fasta subcommand
fasta_example_text = '''
example:
  # cat chrom.csv
  # 1
  # 2
  # 3
  GenoKit fasta -f human.fa -i chrom.csv -o chr1_2_3.fa'''

parser_fasta = subparsers.add_parser('fasta', help='Extract selected FASTA records',
                                     epilog=fasta_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
parser_fasta.add_argument('-f', '--fasta', type=str, required=True,
                           help='fasta file path')
parser_fasta.add_argument('-i', '--input', type=str, required=True,
                           help='input file containing fasta id, one per line')
parser_fasta.add_argument('-o', '--output', type=str, required=True,
                           help='output file path')

parser_fasta.set_defaults(func=fasta)



# IUPAC subcommand https://www.bioinformatics.org/sms/iupac.html
iupac_example_text = '''
example:
  # Amino acid code
  GenoKit iupac -a
  # Codon table
  GenoKit iupac -c
  # Nucleotide code
  GenoKit iupac -n'''
parser_iupac = subparsers.add_parser('iupac', help='Display standard IUPAC codes',
                                     epilog=iupac_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
iupac_options = parser_iupac.add_mutually_exclusive_group()
iupac_options.add_argument('-a', '--amino', action="store_true",
                           help='amino acid code')
iupac_options.add_argument('-c', '--codon', action="store_true",
                           help='codon table')
iupac_options.add_argument('-n', '--nucleotide', action="store_true",
                           help='nucleotide code')
parser_iupac.set_defaults(func=iupac, command_parser=parser_iupac)



# kmer subcommand
kmer_description_text = '''
Build a source-aware k-mer catalogue for GenoKit design commands.

Important k values:
  k=18-25    Complete primer sequences; exact primer lengths are queried.
  k=19-21    Complete sgRNA, siRNA and shRNA target sequences.

Recommended: build k=18-25 for both biological sources. The genome
dataset supports genomic primers and sgRNA; the mature-transcriptome dataset
supports RT-PCR/qPCR primers, siRNA and shRNA. Both sources can coexist in
the same genokit.kmer.sqlite catalogue and are selected independently.
'''
kmer_example_text = '''
example:
  # 1. Genome catalogue: complete primer/sgRNA lengths
  GenoKit kmer -f human.fa --lengths 18-25 -o ./ -r \\
      --source-type genome -j /opt/anaconda3/envs/bio/bin/jellyfish -t 4

  # 2. Add mature-transcript k-mers to the SAME genokit.kmer.sqlite:
  #    complete RT-qPCR/siRNA/shRNA lengths
  GenoKit kmer --source-type transcriptome -d human.gtf.sqlite -g human.fa \\
      --lengths 18-25 -o ./ -r \\
      -j /opt/anaconda3/envs/bio/bin/jellyfish -t 4

  # 3. The resulting catalogue supports, for example:
  GenoKit primer -d human.gtf.sqlite -g human.fa -s gtf -o qpcr.csv \\
      -k genokit.kmer.sqlite --kmer-source transcriptome \\
      --template-type cdna
  GenoKit sgrna -d human.gtf.sqlite -g human.fa -s gtf -o sgrna.csv \\
      -k genokit.kmer.sqlite -l 20
  GenoKit sirna -d human.gtf.sqlite -g human.fa -s gtf -o sirna.csv \\
      -k genokit.kmer.sqlite -l 21
  GenoKit shrna -d human.gtf.sqlite -g human.fa -s gtf -o shrna.csv \\
      -k genokit.kmer.sqlite -l 21

  # 4. Register an existing Jellyfish database without recounting
  GenoKit kmer --register-jf human.fa.kmer.19.jf --source-type genome \\
      -o ./ -r -j /opt/anaconda3/envs/bio/bin/jellyfish
'''
parser_kmer = subparsers.add_parser('kmer', help='Build genome/transcriptome k-mer catalogues',
                                    description=kmer_description_text,
                                    epilog=kmer_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )

# parser_kmer.add_argument('-a', '--annotation', type=str, required=True,
#                           help='genome annotation binary database from create command')
parser_kmer.add_argument('-f', '--fasta', type=str,
                         help='genome FASTA, or a pre-built mature-transcript FASTA')
parser_kmer.add_argument('-d', '--database', type=str,
                         help='annotation database used to build a mature transcriptome')
parser_kmer.add_argument('-g', '--genome', type=str,
                         help='genome FASTA used with --database for transcriptome counting')
parser_kmer.add_argument('--source-type', choices=['genome', 'transcriptome'], default='genome',
                         help='biological source recorded in the k-mer catalogue (default: genome)')
parser_kmer.add_argument('-j', '--jellyfish', type=str,
                          help='path to the jellyfish executable')
parser_kmer.add_argument('--engine', choices=['auto', 'jellyfish', 'python'], default='auto',
                         help='counting engine (default: auto)')
parser_kmer.add_argument('-t', '--threads', type=int,
                         help='worker threads for Jellyfish (default: all available CPUs)')
parser_kmer.add_argument('--hash-size', default='3G',
                         help='Jellyfish hash size (default: 3G)')
parser_kmer.add_argument('-l', '--length', type=int, default=18,
                           help='single k-mer length (default: 18); use --lengths for a design catalogue')
parser_kmer.add_argument('--lengths', type=str,
                         help='multiple k values (recommended: 18-25 for complete primers and RNA targets)')
parser_kmer.add_argument('-o', '--output_dir', type=str, required=True,
                           help='kmer output directory (auto-created)')
parser_kmer.add_argument('--sqlite-db', type=str,
                         help='SQLite k-mer catalogue (default: OUTPUT_DIR/genokit.kmer.sqlite)')
parser_kmer.add_argument('--register-jf', action='append',
                         help='register an existing Jellyfish database; may be repeated')
parser_kmer.add_argument('-r', '--reverse_complement', action="store_true",
                           help='include reverse-complement k-mers (canonical mode with Jellyfish)')

# parser_kmer.add_argument('-p', '--process', type=int, default=4,
#                           help=f'cpu count (default: 4)')
# parser_kmer.add_argument('--single', action='store_true', help='single cpu mode')

parser_kmer.set_defaults(func=kmer)


# stat subcommand
stat_example_text = '''
example:
  # sqlite 
  GenoKit stat -d human.gtf.sqlite -g human.fa -o human.csv -s gtf 
  # pickle
  GenoKit stat -d human.gtf.pkl -g human.fa -o human.csv -s gtf
'''
parser_stat = subparsers.add_parser('stat', help='Summarize annotation statistics',
                                    epilog=stat_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_stat.add_argument('-d', '--database', type=str, required=True,
                           help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_stat.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_stat.add_argument('-o', '--output', type=str,
                           help='output path, if not stdout')
parser_stat.add_argument('-s', '--style', choices=['gff','gtf'],
                           help='genome annotation file format')
parser_stat.set_defaults(func=stat)


# promoter subcommand
promoter_example_text = '''
example:
  # CSV
  GenoKit promoter -d human.gtf.sqlite -f csv -g human.fa -l 1000 -u 10 -o promoter.csv -s gtf
  # GFF 
  GenoKit promoter -d human.gtf.sqlite -f gff -g human.fa -l 1000 -u 10 -o promoter.gff -s gtf
  # GTF
  GenoKit promoter -d human.gtf.sqlite -f gtf -g human.fa -l 1000 -u 10 -o promoter.gtf -s gtf
  # pickle
  GenoKit promoter -d human.gtf.pkl -f csv -g human.fa -l 1000 -u 10 -o promoter.csv -s gtf
  
note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command. 
'''
parser_promoter = subparsers.add_parser('promoter', help='Extract promoter sequence regions',
                                        epilog=promoter_example_text,
                                        formatter_class=argparse.RawTextHelpFormatter
                                        )
parser_promoter.add_argument('-d', '--database', type=str, required=True,
                             help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_promoter.add_argument('-f', '--output_format', type=str,
                             choices=['csv','fasta','gff','gtf'], default='csv',
                             help = 'output format')
parser_promoter.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_promoter.add_argument('-i', '--gene', type=str,
                             help='specific gene id (extract all if omitted)')
parser_promoter.add_argument('-l', '--promoter_length', type=int, default=100,
                             help='promoter length before TSS (default: 100)')
parser_promoter.add_argument('-o', '--output', type=str,
                             help = 'output file path')
# parser_promoter.add_argument('-p', '--process', type=int, default=4,
#                              help='number of promoter extract process, (default: 4)')
parser_promoter.add_argument('-s', '--style', choices=['gff','gtf'],
                             help='genome annotation file format')
parser_promoter.add_argument('-u', '--utr5_upper_length', type=int, default=10,
                             help='5\' utr length after TSS (default: 10)')
parser_promoter.add_argument('-v', '--print', action="store_true",
                             help = 'output to stdout')
parser_promoter.set_defaults(func=promoter)


# terminator
terminator_example_text = '''
example:
  # CSV
  GenoKit terminator -d human.gtf.sqlite -g human.fa -l 1000 -u 10 -o terminator.csv -f csv -s gtf
  # GFF
  GenoKit terminator -d human.gtf.sqlite -g human.fa -l 1000 -u 10 -o terminator.gff -f gff -s gtf
  # GTF
  GenoKit terminator -d human.gtf.sqlite -g human.fa -l 1000 -u 10 -o terminator.gtf -f gtf -s gtf
  # pickle
  GenoKit terminator -d human.gtf.pkl -g human.fa -l 1000 -u 10 -o terminator.csv -f csv -s gtf
  
note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_terminator = subparsers.add_parser('terminator', help='Extract terminator sequence regions',
                                          epilog=terminator_example_text,
                                          formatter_class=argparse.RawTextHelpFormatter
                                          )
parser_terminator.add_argument('-d', '--database', type=str, required=True,
                             help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_terminator.add_argument('-f', '--output_format', type=str, choices=['csv','fasta','gff','gtf'],
                             help = 'output format')
parser_terminator.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_terminator.add_argument('-i', '--gene', type=str,
                             help='specific gene (optional); if not given, return whole genes')
parser_terminator.add_argument('-l', '--terminator_length', type=int, default=100,
                             help='terminator length (default: 100)')
parser_terminator.add_argument('-o', '--output', type=str,
                             help = 'output file path')
parser_terminator.add_argument('-s', '--style', choices=['gff','gtf'],
                             help='genome annotation file format')
parser_terminator.add_argument('-u', '--utr3_lower_length', type=int, default=10,
                             help='3\' length (default: 10)')
parser_terminator.add_argument('-v', '--print', action="store_true",
                             help = 'output to stdout')
parser_terminator.set_defaults(func=terminator)


# gene subcommand
gene_example_text = '''
example:
  # CSV
  GenoKit gene -d human.gtf.sqlite -f csv -g human.fa -o gene.csv -s gtf
  # GFF
  GenoKit gene -d human.gtf.sqlite -f gff -g human.fa -o gene.gff -s gtf
  # Fasta
  GenoKit gene -d human.gtf.sqlite -f fasta -g human.fa -o gene.fasta -s gtf
  # pickle
  GenoKit gene -d human.gtf.pkl -f csv -g human.fa -o gene.csv -s gtf
  GenoKit gene -d human.gtf.pkl -f gtf -g human.fa -o gene.gtf -s gtf
  GenoKit gene -d human.gtf.pkl -f fasta -g human.fa -o gene.fasta -s gtf
  
note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_gene = subparsers.add_parser('gene', help='Extract genomic gene sequences',
                                    epilog=gene_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_gene.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_gene.add_argument('-f', '--output_format', type=str, choices=['csv', 'fasta', 'gff', 'gtf'],
                             help = 'output format')
parser_gene.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_gene.add_argument('-i', '--gene', type=str,
                         help='specific gene id (extract all if omitted)')
parser_gene.add_argument('-o', '--output', type=str,
                         help = 'output file path')
parser_gene.add_argument('-p', '--print', action="store_true",
                         help='output to stdout')
parser_gene.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_gene.set_defaults(func=gene)

# mrna subcommand
# parser_mrna = subparsers.add_parser('mrna', help='extract messager RNA')
# parser_mrna.add_argument('-d', '--database', type=str, required=True,
#                          help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
# parser_mrna.add_argument('-f', '--output_format', type=str, choices=['csv','fasta'], default='csv',
#                          help = 'output format')
# parser_mrna.add_argument('-g', '--genome', type=str, required=True,
#                          help='genome fasta')
# parser_mrna.add_argument('-i', '--transcript', type=str,
#                          help='specific transcript (optional); if not given, return whole transcripts')
# parser_mrna.add_argument('-o', '--output', type=str,
#                          help = 'output file path')
# parser_mrna.add_argument('-p', '--print', action="store_true",
#                          help='output to stdout')
# parser_mrna.add_argument('-s', '--style', choices=['gff','gtf'],
#                          help = 'gtf database or gff database')
# parser_mrna.add_argument('-u', '--upper', action="store_true",
#                          help='upper cds and lower utr')
# parser_mrna.set_defaults(func=mrna)


# transcript subcommand

transcript_example_text = '''
example:
  # CSV 
  GenoKit transcript -d human.gtf.sqlite -f csv -g human.fa -o transcript.csv -s gtf
  # GFF
  GenoKit transcript -d human.gtf.sqlite -f gff -g human.fa -o transcript.gff -s gtf
  # GTF
  GenoKit transcript -d human.gtf.sqlite -f gtf -g human.fa -o transcript.gtf -s gtf
  # CDS and utr
  GenoKit transcript -d human.gtf.sqlite -f fasta -g human.fa -o transcript.fa -s gtf -u
  # pickle
  GenoKit transcript -d human.gtf.pkl -f csv -g human.fa -o transcript.csv -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''

parser_transcript = subparsers.add_parser('transcript', help='Extract spliced transcript sequences',
                                          epilog=transcript_example_text,
                                          formatter_class=argparse.RawTextHelpFormatter
                                          )
parser_transcript.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_transcript.add_argument('-f', '--output_format', type=str, choices=['csv','fasta','gff','gtf'], default='csv',
                         help = 'output format')
parser_transcript.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_transcript.add_argument('-i', '--transcript', type=str,
                         help='specific transcript id (extract all if omitted)')
parser_transcript.add_argument('-o', '--output', type=str,
                         help = 'output file path')
# parser_transcript.add_argument('-p', '--process', type=int, default=4,
#                          help='number of cDNA extract process, (default: 4)')
# parser_transcript.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for extract transcript (default: mrna)')
parser_transcript.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_transcript.add_argument('-u', '--upper', action="store_true",
                         help='upper cds and lower utr')
parser_transcript.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_transcript.set_defaults(func=transcript)


# igr subcommand
igr_example_text = '''
example:
  # CSV
  GenoKit igr -d human.gtf.sqlite -f csv -g human.fa -l 10 -o igr.csv -s gtf
  # GFF
  GenoKit igr -d human.gtf.sqlite -f gff -g human.fa -l 10 -o igr.gff -s gtf
  # GTF 
  GenoKit igr -d human.gtf.sqlite -f gtf -g human.fa -l 10 -o igr.gtf -s gtf 
  # pickle
  GenoKit igr -d human.gtf.pkl -f csv -g human.fa -l 10 -o igr.csv -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_igr = subparsers.add_parser('igr', help='Extract intergenic region sequences',
                                   epilog=igr_example_text,
                                   formatter_class=argparse.RawTextHelpFormatter
                                   )
parser_igr.add_argument('-d', '--database', type=str, required=True,
                        help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_igr.add_argument('-f', '--output_format', type=str, choices=['csv','fasta', 'gff','gtf'], default='csv',
                         help = 'output format')
parser_igr.add_argument('-g', '--genome', type=str, required=True,
                        help='genome fasta')
parser_igr.add_argument('-l', '--igr_length', type=int, default=100,
                        help='igr length threshold')
parser_igr.add_argument('-o', '--output', type=str,
                        help = 'output fasta file path')
# parser_igr.add_argument('-p', '--process', type=int, default=4,
#                          help='number of igr extract process, (default: 4)')
parser_igr.add_argument('-s', '--style', choices=['gff','gtf'],
                        help = 'gtf database or gff database')
parser_igr.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_igr.set_defaults(func=igr)

# utr subcommand
utr_example_text = '''
example:
  # CSV
  GenoKit utr -d human.gtf.sqlite -g human.fa -f csv -o utr.csv -s gtf
  # Fasta
  GenoKit utr -d human.gtf.sqlite -g human.fa -f fasta -o utr.fasta -s gtf
  # GFF
  GenoKit utr -d human.gtf.sqlite -g human.fa -f gff -o utr.gff -s gtf
  # GTF 
  GenoKit utr -d human.gtf.sqlite -g human.fa -f gtf -o utr.gtf -s gtf 
  # pickle 
  GenoKit utr -d human.gtf.pkl -g human.fa -f csv -o utr.csv -s gtf 

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_utr = subparsers.add_parser('utr', help='Extract 5′ and 3′ UTR sequences',
                                   epilog=utr_example_text,
                                   formatter_class=argparse.RawTextHelpFormatter
                                   )
parser_utr.add_argument('-d', '--database', type=str, required=True,
                        help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_utr.add_argument('-f', '--output_format', choices=['csv','fasta','gff','gtf'], default='csv',
                        help='output format (default: csv)')
parser_utr.add_argument('-g', '--genome', type=str, required=True,
                        help='genome fasta file')
parser_utr.add_argument('-i', '--transcript', type=str,
                        help='specific transcript id (extract all if omitted)')
parser_utr.add_argument('-o', '--output', type=str,
                        help='output file path')
# parser_utr.add_argument('-p', '--process', type=int, default=4,
#                          help='number of utr extract process, (default: 4)')
# parser_utr.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for extract utr, (default: mrna)')
parser_utr.add_argument('-s', '--style', choices=['gff','gtf'],
                        help = 'gtf database or gff database')
parser_utr.add_argument('-v', '--print', action="store_true",
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_utr.set_defaults(func=utr)

# uorf subcommand
uorf_example_text = '''
example:
  # CSV
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -f csv -o uorf.csv \\
        -s gtf
        
  # A transcript
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 \\
        -f csv -v -s gtf
        
  # A figure without intron 
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 -f csv \\
        -o ENST00000394705_uorf.csv \\
        -s gtf -m uorf_figure
        
  # A figure with intron 
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 -f csv \\
        -o ENST00000394705_uorf.csv \\
        -s gtf -n uorf_figure_intron
  # GFF
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -f gff -o uorf.gff \\
        -s gtf
  # GTF
  GenoKit uorf -d human.gtf.sqlite -g human.fa \\
        -f gtf -o uorf.gtf \\
        -s gtf -i ENST00000511072
        
  # pickle
  GenoKit uorf -d human.gtf.pkl -g human.fa \\
        -f csv -o uorf.csv \\
        -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''

parser_uorf = subparsers.add_parser('uorf', help='Extract upstream open reading frames',
                                    epilog=uorf_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_uorf.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_uorf.add_argument('-f', '--output_format', choices=['csv','fasta','gff','gtf'], default='csv',
                        help='output format (default: csv)')
parser_uorf.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_uorf.add_argument('-i', '--transcript', type=str,
                         help='specific transcript id (extract all if omitted)')
parser_uorf.add_argument('-l', '--length', type=int, default=6,
                         help='uorf length threshold, (default: 6)')
parser_uorf.add_argument('-m', '--schematic_without_intron', type=str,
                         help=('schematic figure file for uorf, cds\n'
                         'and transcript without intron'))
parser_uorf.add_argument('-n', '--schematic_with_intron', type=str,
                         help=('schematic figure file for uorf, cds\n'
                         'and transcript with intron'))
parser_uorf.add_argument('-o', '--output', type=str,
                         help='output file path')
# parser_uorf.add_argument('-p', '--process', type=int, default=4,
#                          help='number of uorf extract process, (default: 4)')
# parser_uorf.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for uorf extraction (default: mrna)')
parser_uorf.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_uorf.add_argument('-v', '--print', action="store_true",
                         help='output to stdout. -v and -o option are mutually exclusive')
parser_uorf.set_defaults(func=uorf)

# cds subcommand
cds_example_text = '''
example:
  # CSV
  GenoKit cds -d human.gtf.sqlite -f csv -g human.fa -o cds.csv -s gtf
  # Fasta
  GenoKit cds -d human.gtf.sqlite -f fasta -g human.fa -o cds.fasta -s gtf
  # GFF
  GenoKit cds -d human.gtf.sqlite -f gff -g human.fa -o cds.gff -s gtf
  # GTF
  GenoKit cds -d human.gtf.sqlite -f gtf -g human.fa -o cds.gtf -s gtf
  # pickle
  GenoKit cds -d human.gtf.pkl -f csv -g human.fa -o cds.csv -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_cds = subparsers.add_parser('cds', help='Extract coding sequence regions',
                                   epilog=cds_example_text,
                                   formatter_class=argparse.RawTextHelpFormatter
                                   )
parser_cds.add_argument('-d', '--database', type=str, required=True,
                        help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_cds.add_argument('-f', '--output_format', type=str, choices=['csv','fasta','gff','gtf'], default='csv',
                         help = 'output format')
parser_cds.add_argument('-g', '--genome', type=str, required=True,
                        help='genome fasta')
parser_cds.add_argument('-i', '--transcript', type=str,
                        help='specific transcript id (extract all if omitted)')
parser_cds.add_argument('-o', '--output', type=str,
                        help='output file path')

# parser_cds.add_argument('-p', '--process', type=int, default=4,
#                         help='number of cds extract process, (default: 4)')

# parser_cds.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for extract cds, (default: mrna)')
parser_cds.add_argument('-s', '--style', choices=['gff','gtf'],
                        help = 'gtf database or gff database')
parser_cds.add_argument('-v', '--print', action="store_true",
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_cds.set_defaults(func=cds)

# dorf subcommand
dorf_example_text = '''
example:
  # CSV
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -f csv -o dorf.csv \\
        -s gtf
        
  # A transcript
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 \\
        -f csv -v -s gtf
        
  # A figure without intron 
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 -f csv \\
        -o ENST00000394705_dorf.csv \\
        -s gtf -m dorf_figure
        
  # A figure with intron
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -i ENST00000394705 -f csv \\
        -o ENST00000394705_dorf.csv \\
        -s gtf \\
        -n dorf_figure_intron
  # GFF
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -f gff -o dorf.gff \\
        -s gtf
  # GTF
  GenoKit dorf -d human.gtf.sqlite -g human.fa \\
        -f gtf -o dorf.gtf \\
        -s gtf -i ENST00000511072
        
  # pickle
  GenoKit dorf -d human.gtf.pkl -g human.fa \\
        -f csv -o dorf.csv \\
        -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''

parser_dorf = subparsers.add_parser('dorf', help='Extract downstream open reading frames',
                                    epilog=dorf_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_dorf.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_dorf.add_argument('-f', '--output_format', choices=['csv','fasta','gff','gtf'], default='csv',
                         help='output format')
parser_dorf.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_dorf.add_argument('-i', '--transcript', type=str,
                         help='specific transcript id (extract all if omitted)')
parser_dorf.add_argument('-l', '--length', type=int, default=6,
                         help='dorf length, (default: 6)')
parser_dorf.add_argument('-m', '--schematic_without_intron', type=str,
                         help=('schematic figure file for dorf, cds \n'
                         'and transcript without intron'))
parser_dorf.add_argument('-n', '--schematic_with_intron', type=str,
                         help=('schematic figure file for dorf, cds \n'
                         'and transcript with intron'))
parser_dorf.add_argument('-o', '--output', type=str,
                         help='output file path')
# parser_dorf.add_argument('-p', '--process', type=int, default=4,
#                          help='number of dorf extract process, (default: 4)')
# parser_dorf.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for dorf extraction (default: mrna)')
parser_dorf.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_dorf.add_argument('-v', '--print', action="store_true",
                         help='output to stdout. -v and -o option are mutually exclusive')
parser_dorf.set_defaults(func=dorf)


# exon
exon_example_text = '''
example:
  # CSV
  GenoKit exon -d human.gtf.sqlite -g human.fa -f csv -o exon.csv -s gtf
  # GFF
  GenoKit exon -d human.gtf.sqlite -g human.fa -f gff -o exon.gff -s gtf
  # GTF
  GenoKit exon -d human.gtf.sqlite -g human.fa -f gtf -o exon.gtf -s gtf
  # Fasta
  GenoKit exon -d human.gtf.sqlite -g human.fa -f fasta -o exon.fasta -s gtf
  # pickle
  GenoKit exon -d human.gtf.pkl -g human.fa -f csv -o exon.csv -s gtf

note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''

parser_exon = subparsers.add_parser('exon', help='Extract transcript exon sequences',
                                    epilog=exon_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_exon.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_exon.add_argument('-f', '--output_format', choices=['csv','fasta','gff', 'gtf'], default='csv',
                         help='output format')
parser_exon.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_exon.add_argument('-i', '--transcript', type=str,
                         help='specific transcript id (extract all if omitted)')
parser_exon.add_argument('-o', '--output', type=str,
                         help='output file path')
# parser_exon.add_argument('-p', '--process', type=int, default=4,
#                          help='number of exon extract process, (default: 4)')
# parser_exon.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                          help='The type of RNA for exon extraction (default: mrna)')
parser_exon.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_exon.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_exon.set_defaults(func=exon)

# intron
intron_example_text = '''
example:
  # CSV  
  GenoKit intron -d human.gtf.sqlite -g human.fa -f csv -o intron.csv -s gtf
  # GFF
  GenoKit intron -d human.gtf.sqlite -g human.fa -f gff -o intron.gff -s gtf
  # GTF
  GenoKit intron -d human.gtf.sqlite -g human.fa -f gtf -o intron.gtf -s gtf
  # Fasta
  GenoKit intron -d human.gtf.sqlite -g human.fa -f fasta -o intron.fasta -s gtf
  # A transcript
  GenoKit intron -d human.gtf.sqlite -g human.fa -f csv -i ENST00000378486 -o one.csv -s gtf
  
note:
  The human.gtf.sqlite (or human.gtf.pkl) file is generated using the GenoKit create command.
'''
parser_intron = subparsers.add_parser('intron', help='Extract gene intron sequences',
                                      epilog=intron_example_text,
                                      formatter_class=argparse.RawTextHelpFormatter
                                      )
parser_intron.add_argument('-d', '--database', type=str, required=True,
                           help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_intron.add_argument('-f', '--output_format', choices=['csv','fasta','gff','gtf'], default='csv',
                           help='output format')
parser_intron.add_argument('-g', '--genome', type=str, required=True,
                           help='genome fasta')
parser_intron.add_argument('-i', '--transcript', type=str,
                           help='specific transcript id (extract all if omitted)')
parser_intron.add_argument('-o', '--output', type=str,
                           help='output file path')
# parser_intron.add_argument('-p', '--process', type=int, default=4,
#                            help='number of exon extract process, (default: 4)')
# parser_intron.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
#                            help='The type of RNA for intron extraction (default: mrna)')
parser_intron.add_argument('-s', '--style', choices=['gff','gtf'],
                           help = 'gtf database or gff database')
parser_intron.add_argument('-v', '--print', action="store_true",
                           help='output to stdout')
parser_intron.set_defaults(func=intron)

# motif
motif_example_text = '''
example:
  # One pattern
  GenoKit motif -b -c -i example/motif_input.fa -m CTACGAAAG -t dna -o motif.csv
  # Patter file
  # cat pattern.file
  # CTACGAAAG
  GenoKit motif -b -c -i example/motif_input.fa -f pattern.file -t dna -o motif.csv
  # Mismatch 
  GenoKit motif -b -c -i example/motif_input.fa -m CTACGAAAG -t dna -o motif.csv -x 2
  # Regex
  GenoKit motif -b -c -i example/motif_input.fa -m "CTACGAAAG|AAAAAA" -t dna -o motif.csv -x 2 -e
'''

parser_motif = subparsers.add_parser('motif', help='Search and analyze sequence motifs',
                                     epilog=motif_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
parser_motif.add_argument('-b', '--both-strands', action='store_true',
                          help='search on both chains')
parser_motif.add_argument('-c', '--case-sensitive', action='store_true',
                          help='case sensitive')
parser_motif.add_argument('-d', '--degenerate', action='store_true',
                          help='using IUPAC degenerate bases')
parser_motif.add_argument('-e', '--regex', action='store_true',
                          help='using regular expressions')
parser_motif.add_argument('-f', '--file',
                          help='a file containing multiple motifs (one per line)')
parser_motif.add_argument('-i', '--input', type=str, required=True,
                          help='input fasta file')
parser_motif.add_argument('-m', '--motif', type=str,
                          help='sequence pattern')
parser_motif.add_argument('-o', '--output', type=str,
                          help='output file path')
parser_motif.add_argument('-t', '--type', choices=['dna', 'rna'], default='dna',
                          help='seq type: dna or rna')
parser_motif.add_argument('-x', '--mismatches', type=int, default=0,
                          help='Maximum number of mismatches allowed')

# parser_motif.add_argument('-p', '--process', type=int, default=4,
#                            help='number of exon extract process, (default: 4)')
parser_motif.set_defaults(func=motif)

# sirna
sirna_example_text = '''
example:
  GenoKit sirna -d human.gtf.sqlite -f csv -g human.fa -l 21 \\
                -o sirna.csv -s gtf
  GenoKit sirna -d human.gtf.sqlite -f csv -g human.fa -l 21 \\
                -k kmer_db/genokit.kmer.sqlite -o sirna.csv -s gtf
  GenoKit sirna -d yeast.gff.sqlite -f csv -g yeast.fa -l 21 \\
                -k kmer_db/genokit.kmer.sqlite -o sirna.csv -s gff

note:
  • The genokit.kmer.sqlite file is created from the GenoKit kmer command.
  • The human.gtf.sqlite file is generated using the GenoKit create command.
'''

parser_sirna = subparsers.add_parser('sirna', help='Design small interfering RNAs',
                                     epilog=sirna_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
parser_sirna.add_argument('-d', '--database', type=str, required=True,
                          help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_sirna.add_argument('-f', '--output_format', choices=['csv','fasta'], default='csv',
                          help='output format (default: csv)')
parser_sirna.add_argument('-g', '--genome', type=str, required=True,
                          help='genome fasta file')
parser_sirna.add_argument('-i', '--gene', type=str,
                          help='specific gene id (extract all if omitted)')
parser_sirna.add_argument('-k', '--kmer_bin', type=str,
                           help='source-aware GenoKit k-mer SQLite database')
parser_sirna.add_argument('--jellyfish', type=str,
                          help='optional Jellyfish executable for .jf queries; specify if the saved path is invalid or Jellyfish is not on PATH')
parser_sirna.add_argument('-l', '--length', type=int, choices=[19, 20, 21], default=21,
                          help='siRNA length (default: 21)')
parser_sirna.add_argument('--region', choices=['cds', 'transcript'], default='cds',
                          help='target region: CDS (default, compatible) or mature transcript')
parser_sirna.add_argument('--target-mode', choices=['all', 'any'], default='all',
                          help='target all input transcripts of each gene (default) or at least one transcript')
parser_sirna.add_argument('--top', type=int,
                          help='maximum number of ranked siRNAs to keep per gene')
parser_sirna.add_argument('-o', '--output', type=str,
                          help='output file path')
# parser_sirna.add_argument('-p', '--process', type=int, default=4,
#                           help='number of sirna process, (default: 4)')
parser_sirna.add_argument('-s', '--style', choices=['gff','gtf'],
                          help = 'gtf database or gff database')
parser_sirna.add_argument('-v', '--print', action="store_true",
                          help='output to stdout. -v and -o option are mutually exclusive')
parser_sirna.set_defaults(func=sirna)

# primer
primer_example_text = '''
example:
  GenoKit primer -a 0.4 -b 0.6  \\
        -c 50 -e 70 -d human.gtf.sqlite \\
        -f csv -g human.fa -m 25 -n 18 \\
        -o primer.csv -s gtf \\
        -x 100 -y 1000 \\
        --template-type genomic --kmer-source genome \\
        -k kmer_db/genokit.kmer.sqlite -l 19

  # Genomic DNA template (build genome primer-length datasets first)
  GenoKit primer -d human.gtf.sqlite -g human.fa -o primer.csv -s gtf \\
      --kmer-source genome \\
      -k kmer_db/genokit.kmer.sqlite \\
      --template-type genomic

  # Mature cDNA template uses transcriptome-labelled k-mer datasets
  GenoKit primer -d human.gtf.sqlite -g human.fa -o qpcr.csv -s gtf \\
      --kmer-source transcriptome \\
      -k kmer_db/genokit.kmer.sqlite \\
      --template-type cdna

note:
  • The genokit.kmer.sqlite is created from the GenoKit kmer command.
  • The human.gtf.sqlite file is generated using the GenoKit create command.
'''
parser_primer = subparsers.add_parser('primer', help='Design gene-specific primer pairs',
                                      epilog=primer_example_text,
                                      formatter_class=argparse.RawTextHelpFormatter
                                      )
parser_primer.add_argument('-a', '--min_gc', type=float, default=0.4,
                           help='min GC content, (default: 0.4)')
parser_primer.add_argument('-b', '--max_gc', type=float, default=0.6,
                           help='max GC content, (default: 0.6)')
parser_primer.add_argument('-c', '--min_tm', type=int, default=50,
                           help='min tm value, (default: 50)')
parser_primer.add_argument('-e', '--max_tm', type=int, default=70,
                           help='max tm value, (default: 70)')
parser_primer.add_argument('-d', '--database', type=str, required=True,
                           help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_primer.add_argument('-f', '--output_format', choices=['fasta','csv'], default='csv',
                           help='output format (default: csv)')
parser_primer.add_argument('-g', '--genome', type=str, required=True,
                           help='genome fasta file')
parser_primer.add_argument('-i', '--gene', type=str,
                           help='specific gene (optional); if not given,\n'
                           'return all genes in database')

parser_primer.add_argument('-k', '--kmer_bin', type=str,
                           help='GenoKit k-mer SQLite database used to calculate Kmer_freq')

parser_primer.add_argument('-l', '--kmer_size', type=int,
                           help='the kmer length of kmer binary database ')
parser_primer.add_argument('--kmer-source', choices=['genome', 'transcriptome'],
                           help='k-mer source; defaults to genome for genomic templates or transcriptome for cdna; must match --template-type')

parser_primer.add_argument('-m', '--max_length', type=int, default=25,
                           help='primer length (bp), (default: 25)')
parser_primer.add_argument('-n', '--min_length', type=int, default=18,
                           help='primer length (bp), (default: 18)')

parser_primer.add_argument('-o', '--output', type=str,
                           help='output file path')

# parser_primer.add_argument('-p', '--process', type=int, default=4,
#                          help='number of primer process, (default: 4)')

parser_primer.add_argument('-s', '--style', choices=['gff','gtf'],
                           help = 'gtf database or gff database')

parser_primer.add_argument('-x', '--min_size', type=int, default=100,
                           help='min target region length (bp), (default: 100)')
parser_primer.add_argument('-y', '--max_size', type=int, default=1000,
                           help='max target region length (bp), (default: 1000)')

parser_primer.add_argument('-v', '--print', action="store_true",
                           help='output to stdout. -v and -o option are mutually exclusive')
parser_primer.add_argument('--template-type', choices=['genomic', 'cdna'], default='cdna',
                           help='primer template: genomic gene DNA including introns, or mature spliced cdna (default: cdna)')
parser_primer.add_argument(
    '--jellyfish',
    help='optional Jellyfish executable for .jf queries; specify if the saved path is invalid or Jellyfish is not on PATH',
)

parser_primer.set_defaults(func=primer)

# sgrna
sgrna_example_text = '''
example:
  # Show PAM table
  GenoKit sgrna --pam
  
  # sgRNA
  GenoKit sgrna -d human.gtf.sqlite -g human.fa -f csv -p NGG -s gtf \\
                -k kmer_db/genokit.kmer.sqlite -o sgrna.csv
  # Kmer
  GenoKit sgrna -d human.gtf.sqlite -g human.fa -f csv -p NGG -s gtf \\
                -k kmer_db/genokit.kmer.sqlite -o sgrna.csv
  # Jellyfish path is missing in SQLite
  GenoKit sgrna -d yeast.gff.sqlite -g yeast.fa -f csv -p NGG -s gff \\
                -k kmer_db/genokit.kmer.sqlite -o sgrna.csv \\
                --jellyfish /opt/anaconda3/envs/bio/bin/jellyfish
note:
  • The genokit.kmer.sqlite file is created from the GenoKit kmer command.
  • The human.gtf.sqlite file is generated using the GenoKit create command.
'''
parser_sgrna = subparsers.add_parser('sgrna', help='Design CRISPR single-guide RNAs',
                                     epilog=sgrna_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
parser_sgrna.add_argument('-d', '--database', type=str,
                          # required=True,
                          help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_sgrna.add_argument('-f', '--output_format', choices=['csv','fasta'], default='csv',
                          help='output format (default: csv)')
parser_sgrna.add_argument('-g', '--genome', type=str,
                          # required=True,
                          help='genome fasta file')
parser_sgrna.add_argument('-i', '--gene', type=str,
                          help='specific gene id (extract all if omitted)')
parser_sgrna.add_argument('-k', '--kmer_bin', type=str,
                           help='genome-source GenoKit k-mer SQLite database')
parser_sgrna.add_argument('--jellyfish', type=str,
                          help='optional Jellyfish executable for .jf queries; specify if the saved path is invalid or Jellyfish is not on PATH')
parser_sgrna.add_argument('-l', '--length', type=int, default=20,
                          help='sgRNA length (default: 20)')
parser_sgrna.add_argument('-o', '--output', type=str,
                          help='output file path')
parser_sgrna.add_argument('-p', '--pam_pattern', type=str, default='NGG',
                          help='PAM pattern, default: NGG (S.pyogenes Cas9)')
parser_sgrna.add_argument('-s', '--style', choices=['gff','gtf'],
                          help = 'gtf database or gff database')
parser_sgrna.add_argument('-v', '--print', action="store_true",
                          help='output to stdout. -v and -o option are mutually exclusive')
parser_sgrna.add_argument('--pam', action="store_true",
                          help='PAM (protospacer adjacent motif) table')

parser_sgrna.set_defaults(func=sgrna)


# shrna
shrna_example_text = '''
example:
  GenoKit shrna -d human.gtf.sqlite -g human.fa -f csv -p CTCGAG -s gtf \\
                -o shrna.csv
  GenoKit shrna -d human.gtf.sqlite -g human.fa -f csv -p CTCGAG -s gtf \\
                -k kmer_db/genokit.kmer.sqlite -o shrna.csv
  GenoKit shrna -d yeast.gff.sqlite -g yeast.fa -f csv -p CTCGAG -s gff \\
                -k kmer_db/genokit.kmer.sqlite -o shrna.csv
                
note:
  • The genokit.kmer.sqlite file is created from the GenoKit kmer command.
  • The human.gtf.sqlite file is generated using the GenoKit create command.
'''

parser_shrna = subparsers.add_parser('shrna', help='Design small hairpin RNA candidates',
                                     epilog=shrna_example_text,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
parser_shrna.add_argument('-d', '--database', type=str, required=True,
                          help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_shrna.add_argument('-f', '--output_format', choices=['csv','fasta'], default='csv',
                          help='output format (default: csv)')
parser_shrna.add_argument('-g', '--genome', type=str, required=True,
                          help='genome fasta file')
parser_shrna.add_argument('-i', '--gene', type=str,
                          help='specific gene id (extract all if omitted)')
parser_shrna.add_argument('-k', '--kmer_bin', type=str,
                           help='transcriptome-source GenoKit k-mer SQLite database')
parser_shrna.add_argument('--jellyfish', type=str,
                          help='optional Jellyfish executable for .jf queries; specify if the saved path is invalid or Jellyfish is not on PATH')
parser_shrna.add_argument('-l', '--length', type=int, choices=[19, 20, 21], default=21,
                          help='shRNA target length: 19, 20, or 21 nt (default: 21)')
parser_shrna.add_argument('-o', '--output', type=str,
                          help='output file path')
parser_shrna.add_argument('-p', '--loop', type=str, default='CTCGAG',
                          help='loop: an ATCG sequence, or standard/miR30/simple (default: CTCGAG)')
parser_shrna.add_argument('-s', '--style', choices=['gff','gtf'],
                          help = 'gtf database or gff database')
parser_shrna.add_argument('--region', choices=['cds', 'transcript'], default='cds',
                          help='target region: CDS (default; exons if no CDS) or mature transcript')
parser_shrna.add_argument('--target-mode', choices=['all', 'any'], default='all',
                          help='target all input transcripts of each gene (default) or at least one transcript')
parser_shrna.add_argument('-v', '--print', action="store_true",
                          help='output to stdout. -v and -o option are mutually exclusive')

parser_shrna.set_defaults(func=shrna)

# barcode
barcode_example_text = '''
example:
  GenoKit barcode -f csv -l 8 -s dna -o barcode_dna.csv
  GenoKit barcode -f csv -l 8 -s rna -o barcode_rna.csv
  GenoKit barcode -f csv -l 8 -s dna -o barcode_dna.fasta
  GenoKit barcode -f csv -l 8 -s rna -o barcode_rna.fasta
'''
parser_barcode = subparsers.add_parser('barcode', help='Design DNA and RNA barcode sets',
                                       epilog=barcode_example_text,
                                       formatter_class=argparse.RawTextHelpFormatter
                                       )
parser_barcode.add_argument('-f', '--output_format', choices=['csv','fasta'], default='csv',
                          help='output format (default: csv)')


parser_barcode.add_argument('-l', '--length', type=int, default=8, required=True,
                          help='barcode length (default: 8)')

parser_barcode.add_argument('-n', '--count', type=int, default=20,
                          help='the number of barcode (default: 20)')

parser_barcode.add_argument('-o', '--output', type=str,
                          help='output file path')

parser_barcode.add_argument('-s', '--style', choices=['dna','rna'],
                          help = 'style of nucleotide (dna, rna)')

parser_barcode.add_argument('-v', '--print', action="store_true",
                          help='output to stdout. -v and -o option are mutually exclusive')

parser_barcode.set_defaults(func=barcode)

# view
view_example_text = '''
example:
  GenoKit view -i ENSG00000138381 \\
               -d human.gtf.sqlite \\
               -l 6 -f pdf -g human.fa \\
               -o ./ -s gtf --orf-label-limit 20
  
note:
  • The human.gtf.sqlite file is generated using the GenoKit create command.
'''
parser_view = subparsers.add_parser('view', help='Visualize annotated gene structure',
                                    epilog=view_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_view.add_argument('-d', '--database', type=str, required=True,
                           help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_view.add_argument('-f', '--output_format', choices=['pdf','png','tiff'], default='pdf',
                           help='output format (default: pdf)')
parser_view.add_argument('-g', '--genome', type=str, required=True,
                           help='genome fasta file')
parser_view.add_argument('-i', '--gene', type=str,
                           help='specific gene (optional); if not given, return all genes in database')
parser_view.add_argument('-l', '--length', type=int, default=6,
                           help='minimum uORF/dORF length in nt; shorter ORFs are hidden (default: 6)')
parser_view.add_argument('--orf-label-limit', type=int, default=20,
                           help=('show individual ORF numbers only when a transcript has no more '
                                 'than this many ORFs; use 0 to hide all numbers (default: 20)'))
parser_view.add_argument('-o', '--output_dir', type=str,
                           help='output directory path')
# parser_view.add_argument('-p', '--process', type=int, default=4,
#                            help='number of sirna process, (default: 4)')
parser_view.add_argument('-s', '--style', choices=['gff','gtf'],
                           help = 'gtf database or gff database')
parser_view.add_argument('-v', '--print', action="store_true",
                           help='output to stdout. -v and -o option are mutually exclusive')
parser_view.set_defaults(func=view)

# circos
circos_example_text = '''
example:
  # Basic usage with multiple tracks
  GenoKit circos \\
      -a example/chromosome_length.csv \\
      -b example/barplot_data.csv \\
      -c example/cytoband_data.csv \\
      -k example/link_data.csv \\
      -m example/heatmap_data.csv \\
      -p example/point_data.csv \\
      -o ./circos_plot/global.pdf \\
      -f pdf
      
  # Minimal usage (only chromosome and barplot)
  GenoKit circos \\
      -a example/chromosome_length.csv \\
      -b example/barplot_data.csv \\
      -o ./circos_plot/cb.pdf
  
  # Three tracks
  GenoKit circos \\
      -a example/chromosome_length.csv \\ 
      -b example/barplot_data.csv \\
      -k example/link_data.csv \\
      -o ./circos_plot/cbl.pdf

  # Output to stdout (no image files created)
  GenoKit circos \\
      -a example/chromosome_length.csv \\
      -b example/barplot_data.csv \\
      -v -f png > cp.png

  Notes:
  • All input files should be comma-separated (CSV) or tab-separated.
  • Example coordinates use 1-based inclusive start and end positions.
  • Missing relative example/... paths are resolved from GenoKit's installed example data.
  • The -t option automatically distinguishes line and link tracks from the header.
  • Selected tracks are packed automatically from outside to inside: point, line, bar, heatmap.
  • The -v (print) option writes the plot to stdout and is exclusive with -o.
  • Use -f to specify output format: pdf (default), png, or tiff.
'''
parser_circos = subparsers.add_parser('circos', help='Visualize circular genome structure',
                                      epilog=circos_example_text,
                                      formatter_class=argparse.RawTextHelpFormatter)
parser_circos.add_argument('-a', '--chrom', type=str,
                           help="chromosome length file (chr,start,end)")
parser_circos.add_argument('-b', '--bar', type=str,
                           help="barplot track file (chr,start,end,value)")
parser_circos.add_argument('-c', '--cytoband', type=str,
                           help="cytoBand track file (chr1,1,2300000,p36.33,gneg)")
# parser_circos.add_argument('-d', '--database', type=str, required=True,
#                            help='SQLite (.sqlite) or pickle (.pkl) database generated by create')
parser_circos.add_argument('-f', '--output_format', choices=['pdf','png','tiff'], default='pdf',
                           help='output format (default: pdf)')
# parser_circos.add_argument('-g', '--genome', type=str, required=True,
#                            help='genome fasta file')
parser_circos.add_argument('-k', '--link', type=str,
                           help='link track file (chr1,start1,end1,chr2,start2,end2)')
parser_circos.add_argument('-l', '--line', type=str,
                           help='line track file (chr,start,end,value)')
parser_circos.add_argument('-t', '--track', type=str,
                           help='line or link track file (automatically detected from its header)')
parser_circos.add_argument('-m', '--heatmap', type=str,
                           help='heatmap track file (chr,start,end,value)')
parser_circos.add_argument('-w', '--window', type=int, default=1000000,
                           help='window size, (default: 1000000)')
circos_output = parser_circos.add_mutually_exclusive_group(required=True)
circos_output.add_argument('-o', '--output_file', type=str,
                           help='output file path')
parser_circos.add_argument('-p', '--point', type=str,
                            help='Point track file (chr,start,end,value)')
parser_circos.add_argument('-s', '--style', choices=['gff','gtf'],
                           help = 'gtf database or gff database')
circos_output.add_argument('-v', '--print', action="store_true",
                           help='output to stdout. -v and -o option are mutually exclusive')
parser_circos.set_defaults(func=circos)



ppi_example_text = '''
example:
  GenoKit ppi -i example/ppi_test.csv -s tab -o ppi_dir/ --layout forceatlas2
  GenoKit ppi -i example/ppi_test.csv -s tab -o ppi_dir/ --layout star
  GenoKit ppi -i example/ppi_test.csv -s tab -o ppi_dir/ --layout fr
'''
parser_ppi = subparsers.add_parser('ppi', help='Visualize protein interaction networks',
                                   epilog=ppi_example_text,
                                   formatter_class=argparse.RawTextHelpFormatter
                                   )
parser_ppi.add_argument('-i', '--input', type=str, required=True,
                           help='PPI input file (proteinA,proteinB)')
parser_ppi.add_argument('-s', '--separator', choices=['tab','comma'],
                           help='separator between proteinA and proteinB')
parser_ppi.add_argument('-o', '--output_dir', type=str,
                           help='output directory path')
parser_ppi.add_argument('-l', '--layout',
                           choices=['auto', 'fr', 'kk', 'forceatlas2', 'circle',
                                    'star', 'tree', 'grid', 'shell', 'spectral',
                                    'spiral', 'random'],
                           default='auto',
                           help=('network layout inspired by igraph; auto chooses a layout '
                                 'from network topology and size (default: auto)'))
parser_ppi.set_defaults(func=ppi)


hgvs_description_text = '''
Build a tool for visualizing HGVS nomenclature.
A letter prefix indicates the type of reference sequence used.
Accepted prefixes are;
    • c for a coding DNA reference sequence,
    • g for a linear genomic reference sequence,
    • m for a mitochondrial DNA reference sequence,
    • n for a non-coding DNA reference sequence,
    • o for a circular genomic reference sequence,
    • p for a protein reference sequence,
    • r for an RNA reference sequence (transcript).
'''

hgvs_example_text = '''
examples:
  # sub: substitution (DNA, RNA, protein)
  GenoKit hgvs -s "NM_004006.3:c.123C>G" -o dna.png
  GenoKit hgvs -s "NM_004006.3:r.123c>g" -o rna.png
  GenoKit hgvs -s "NP_003997.1:p.Trp24Cys" -o protein.png

  # del: deletion
  GenoKit hgvs -s "NC_000001.11:g.1234_2345del" -o del_dna.png
  GenoKit hgvs -s "NM_004006.3:r.123_127del" -o indel_rna.png

  # ins: insertion (the two positions must be adjacent)
  GenoKit hgvs -s "NC_000001.11:g.1234_1235insACGT" -o ins_dna.png
  GenoKit hgvs -s "NM_004006.3:r.123_124insg" -o ins_rna.png
  GenoKit hgvs -s "NM_004006.3:r.123_124insauc" --skip-three-prime-check -o ins_rna.png

  # delins: deletion-insertion
  GenoKit hgvs -s "NC_000001.11:g.123_129delinsAC" -o delins_dna.png
  GenoKit hgvs -s "NP_004371.2:p.(Asn47delinsSerSerTer)" -o delins_protein.png

  # dup: tandem duplication
  GenoKit hgvs -s "NM_004006.3:r.124_346dup" -o dup_rna.png --skip-three-prime-check
  GenoKit hgvs -s "NP_003997.1:p.Lys23_Val25dup" -o dup_protein.png

  # inv: inversion (DNA/RNA only; at least two nucleotides)
  GenoKit hgvs -s "NC_000001.11:g.1234_2345inv" -o inv_dna.png
  GenoKit hgvs -s "NM_004006.3:r.123_345inv" -o inv_rna.png

  # rpt: repeated sequence
  GenoKit hgvs -s "NC_000014.8:g.19000079CAG[23]" -o repeat_dna.png
  GenoKit hgvs -s "NM_004006.3:r.9495_9497[4]" -o repeat_rna.png
  GenoKit hgvs -s "NP_003997.1:p.Leu2[10]" -o repeat_protein.png

  # alleles: cis, trans, or unknown phase
  GenoKit hgvs -s "NC_000001.11:g.[123G>A;345del]" -o cis.png --allow-reference-mismatch
  GenoKit hgvs -s "NC_000001.11:g.[123G>A];[345del]" -o trans.png --allow-reference-mismatch
  GenoKit hgvs -s "NM_004006.3:r.123c>a(;)345del" -o phase_unknown.png --skip-three-prime-check

  # ext and fs: protein extension and frameshift
  GenoKit hgvs -s "NP_003997.1:p.Met1ext-5" -o extension.png
  GenoKit hgvs -s "NP_003997.1:p.Thr97ProfsTer23" -o frameshift.png
  
note:
  • numbering:https://hgvs-nomenclature.org/stable/background/numbering/
'''
parser_hgvs = subparsers.add_parser('hgvs', help='Visualize DNA, RNA, and protein variants',
                                    description=hgvs_description_text,
                                    epilog=hgvs_example_text,
                                    formatter_class=argparse.RawTextHelpFormatter
                                    )
parser_hgvs.add_argument('-i', '--input', type=str,
                           help='text/CSV file containing one HGVS expression per line')
parser_hgvs.add_argument('-s', '--string',
                           help='DNA, RNA or protein HGVS expression (g./c./n./m./o., r. or p.)')
parser_hgvs.add_argument('--variation',
                           choices=['sub', 'del', 'ins', 'delins', 'dup', 'inv',
                                    'rpt', 'alleles', 'ext', 'fs'],
                           help='optional expected variation type used for validation')
parser_hgvs.add_argument('-o', '--output', type=str,
                           help='output PNG/PDF/SVG path; for multiple inputs, output directory or prefix')
parser_hgvs.add_argument('--genbank', type=str,
                           help='local GenBank/GenPept file instead of downloading from NCBI')
parser_hgvs.add_argument('--cache-dir', default='.genokit_cache/hgvs',
                           help='NCBI GenBank cache directory (default: .genokit_cache/hgvs)')
parser_hgvs.add_argument('--context', type=int, default=20,
                           help='flanking bases/residues displayed in the detailed panel (default: 20)')
parser_hgvs.add_argument('--dpi', type=int, default=200,
                           help='raster output resolution (default: 200)')
parser_hgvs.add_argument('--email', type=str,
                           help='email sent to NCBI E-utilities (or set NCBI_EMAIL)')
parser_hgvs.add_argument('--api-key', type=str,
                           help='NCBI API key (or set NCBI_API_KEY)')
parser_hgvs.add_argument('--allow-reference-mismatch', action='store_true',
                           help='draw an invalid reference allele as a warning instead of stopping')
parser_hgvs.add_argument('--skip-three-prime-check', action='store_true',
                           help="skip the local 3' check for simple del/dup/ins diagnostic input")
parser_hgvs.set_defaults(func=hgvs)

# Reserved for the future expression-pattern command.  It deliberately has no
# dispatcher entry or handler until the implementation is complete.
parser_expp = subparsers.add_parser('expp', help='Expression pattern')
parser_expp.add_argument('-d', '--database', type=str, required=True,
                         help='SQLite (.sqlite) or pickle (.pkl) database generated by create')

# With no options, argparse would validate required primer options before the
# command handler can run.  Treat the bare primer command as a help request.
if len(sys.argv) == 2 and sys.argv[1] == 'primer':
    parser_primer.print_help()
    raise SystemExit(0)

_BARE_COMMAND_HELP_PARSERS = {
    'create': parser_create,
    'fasta': parser_fasta,
    'stat': parser_stat,
    'kmer': parser_kmer,
    'iupac': parser_iupac,
    'gene': parser_gene,
    'transcript': parser_transcript,
    'igr': parser_igr,
    'promoter': parser_promoter,
    'sirna': parser_sirna,
    'sgrna': parser_sgrna,
    'primer': parser_primer,
    'motif': parser_motif,
    'ppi': parser_ppi,
    'view': parser_view,
    'circos': parser_circos,
    'terminator': parser_terminator,
    'utr': parser_utr,
    'uorf': parser_uorf,
    'cds': parser_cds,
    'dorf': parser_dorf,
    'exon': parser_exon,
    'intron': parser_intron,
    'barcode': parser_barcode,
    'shrna': parser_shrna,
    'hgvs': parser_hgvs,
}

if len(sys.argv) == 2 and sys.argv[1] in _BARE_COMMAND_HELP_PARSERS:
    _BARE_COMMAND_HELP_PARSERS[sys.argv[1]].print_help()
    raise SystemExit(0)

args = parser.parse_args()

#print('[%s runing ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))
args.func(args)
#print('[%s finished ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))
