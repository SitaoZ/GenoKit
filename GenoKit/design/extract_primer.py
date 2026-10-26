# -*- coding: utf-8 -*-
"""Design primers on genomic DNA or mature cDNA templates.

CSV/FASTA Start and End are 1-based inclusive coordinates on the oriented
template, for both primer directions. Internal design coordinates stay
0-based half-open; these are not absolute genomic coordinates.
"""

import csv
import logging
import sys
from contextlib import ExitStack

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import get_gene_symbol, load_annotation_database
from GenoKit.database.kmer_database import open_kmer_database
from GenoKit.design.design_primer import PrimerDesigner

# 1-based coordinate [Start, End]
_CSV_HEADER = [
    'Gene ID', 'Gene Symbol', 'Transcript ID', 'Template Type', 'Source',
    'Primer', 'Start', 'End', 'GC', 'Tm', 'Length', 'Fragment',
    'Kmer Freq', 'Template Seq'
]
_TEMPLATE_SOURCES = {'genomic': 'genome', 'cdna': 'transcriptome'}

# Shared ordering keeps CSV columns and FASTA metadata consistent. Templates
# mean all isoforms of the target gene for cDNA, or its one genomic interval.
# Target Template Count: 目标模板数量
# Target Templates Hit:命中数量及Target Coverage:覆盖率
# Target Exact Matches: 目标内精确匹配次数
# Off-target Exact Matches: 非目标精确匹配次数及Specificity Status: 计数状态
# Target Pair Templates Hit: 同一对引物的目标模板覆盖数量及Target Pair Coverage:覆盖率
_SPECIFICITY_FIELDS = (
    ('Target Template Count', 'target_template_count'),
    ('Target Templates Hit', 'target_templates_hit'),
    ('Target Coverage', 'target_coverage'),
    ('Target Exact Matches', 'target_exact_matches'),
    ('Off-target Exact Matches', 'off_target_exact_matches'),
    ('Specificity Status', 'specificity_status'),
    ('Target Pair Templates Hit', 'target_pair_templates_hit'),
    ('Target Pair Coverage', 'target_pair_coverage'),
)
_CSV_HEADER.extend(label for label, _ in _SPECIFICITY_FIELDS)


def _specificity_values(primer, pair):
    """Coverage is a fraction in [0, 1]; unavailable counts remain explicit NA."""
    return [primer.get(key, pair.get(key, 'NA')) for _, key in _SPECIFICITY_FIELDS]


def _primer_cli_error(message):
    raise SystemExit(f"GenoKit primer: error: {message}")


def parse_transcript(params):
    """Return strand-oriented, spliced templates keyed by transcript ID."""
    g, ts, db, genome, style, output_format = params
    templates = {}
    for t in ts:
        model = db[g][t]
        sequence = ''.join(
            str(genome[e.chr][e.start - 1:e.end])
            for e in sorted(model['exon'], key=lambda e: (e.start, e.end))
        )
        if model['mRNA'].strand == '-':
            sequence = str(Seq(sequence).reverse_complement())
        templates[t] = sequence.upper()
    return {g: templates}


def primer_templates(gene_id, db, genome, template_type):
    """Genomic DNA uses one continuous gene interval; cDNA uses each isoform.

    Both templates are oriented in the gene/transcript's 5-prime direction.
    Primer coordinates are zero-based, half-open offsets into these templates.
    """
    gene = db[gene_id]
    if template_type == 'genomic':
        feature = gene['gene']
        sequence = str(genome[feature.chr][feature.start - 1:feature.end])
        if feature.strand == '-':
            sequence = str(Seq(sequence).reverse_complement())
        return {gene_id: sequence.upper()}
    # cdna
    ts = [t for t in gene if t not in {'gene', 'symbol'}]
    return parse_transcript((gene_id, ts, db, genome, None, None))[gene_id]


def _template_source(args):
    expected = _TEMPLATE_SOURCES[args.template_type]
    if args.kmer_source and args.kmer_source != expected:
        _primer_cli_error(
            f"--template-type {args.template_type} requires --kmer-source {expected}; "
            f"received {args.kmer_source}."
        )
    if args.jellyfish and not args.kmer_bin:
        _primer_cli_error('--jellyfish requires a k-mer database supplied with -k/--kmer_bin.')
    return expected


def _close_when_done(stack, resource):
    close = getattr(resource, 'close', None)
    if close is not None:
        stack.callback(close)
    return resource


def get_primer(args):
    expected_source = _template_source(args)
    if not args.print and not args.output:
        _primer_cli_error('provide -o/--output or --print.')
    if args.print and args.output:
        _primer_cli_error('use either -o/--output or --print, not both.')
    if not (0 <= args.min_gc <= args.max_gc <= 1):
        _primer_cli_error('GC bounds must satisfy 0 <= min_gc <= max_gc <= 1.')
    if not (0 < args.min_length <= args.max_length):
        _primer_cli_error('primer lengths must satisfy 0 < min_length <= max_length.')
    if not (0 < args.min_size <= args.max_size):
        _primer_cli_error('product sizes must satisfy 0 < min_size <= max_size.')
    if args.min_tm > args.max_tm:
        _primer_cli_error('min_tm must not exceed max_tm.')

    logging.basicConfig(level=logging.INFO, filename='Primer.log')
    with ExitStack() as stack:
        try:
            counts = open_kmer_database(
                args.kmer_bin, k=args.kmer_size,
                expected_source=expected_source, jellyfish=args.jellyfish,
            )
        except (FileNotFoundError, ValueError) as error:
            _primer_cli_error(str(error))
        _close_when_done(stack, counts)
        available = tuple(getattr(counts, 'available_k', ()))
        if args.kmer_bin and available:
            lengths = set(range(args.min_length, args.max_length + 1))
            supported = lengths.intersection(available)
            if not supported:
                _primer_cli_error(
                    f'The {expected_source} k-mer datasets contain k={available}, '
                    f'but primer lengths are {args.min_length}-{args.max_length}. '
                    'Build matching exact k values with GenoKit kmer --lengths.'
                )
            missing = sorted(lengths.difference(available))
            if missing:
                print(f'Warning: missing exact k values {missing}; primer lengths '
                      f'are restricted to {sorted(supported)}.', file=sys.stderr)

        db, t2g = load_annotation_database(args.database)
        _close_when_done(stack, db)
        _close_when_done(stack, t2g)
        genome = _close_when_done(stack, genome_dict(args.genome))
        if args.gene and args.gene not in db:
            _primer_cli_error(f'gene ID {args.gene} does not exist in the annotation database.')
        output = sys.stdout if args.print else stack.enter_context(
            open(args.output, 'w', newline=''))
        writer = None
        if args.output_format == 'csv':
            writer = csv.writer(output)
            writer.writerow(_CSV_HEADER)
        genes = [args.gene] if args.gene else db
        for g in tqdm(genes, total=1 if args.gene else len(db),
                      desc='Primer processing:', ncols=80, colour="GREEN"):
            templates = primer_templates(g, db, genome, args.template_type)
            designer = PrimerDesigner(
                g, templates, counts, args.kmer_size or 0,
                min_gc=args.min_gc, max_gc=args.max_gc,
                min_len=args.min_length, max_len=args.max_length,
                min_tm=args.min_tm, max_tm=args.max_tm,
                min_size=args.min_size, max_size=args.max_size,
                template_type=args.template_type,
            )
            for template_id in templates:
                designer.design_primers_for_transcript(template_id)
            primers = designer.get_all_primer_info()
            if not primers:
                if args.gene:
                    print(f'Warning: no primer pair met the requested constraints for {g}.',
                          file=sys.stderr)
                continue
            if args.gene:
                designer.visualize_primers_multi(list(primers))
            symbol = get_gene_symbol(db, g)
            for template_id, pair in primers.items():
                transcript_id = template_id if args.template_type == 'cdna' else ''
                for direction in ('forward', 'reverse'):
                    primer = pair[direction]
                    # Convert only at serialization: [start, end) -> [start+1, end].
                    # The reverse primer also uses ascending template coordinates.
                    output_start, output_end = primer['start'] + 1, primer['end']
                    if writer is not None:
                        writer.writerow([
                            g, symbol, transcript_id, args.template_type, direction,
                            primer['sequence'], output_start, output_end,
                            round(primer['gc'], 2), primer['tm'], primer['length'],
                            pair['product_size'], primer.get('kmer_freq', 'NA'),
                            templates[template_id],
                        ] + _specificity_values(primer, pair))
                    else:
                        specificity = ' '.join(
                            f'{key}:{value}' for (_, key), value in zip(
                                _SPECIFICITY_FIELDS, _specificity_values(primer, pair)))
                        record = SeqRecord(
                            Seq(primer['sequence']), id=f'{g}|{template_id}|{direction}',
                            description=(f'template_type:{args.template_type} '
                                         f'start:{output_start} end:{output_end} '
                                         f'length:{primer["length"]} '
                                         f'product_size:{pair["product_size"]} '
                                         f'kmer_freq:{primer.get("kmer_freq", "NA")} '
                                         f'{specificity}'),
                        )
                        SeqIO.write(record, output, 'fasta')


if __name__ == "__main__":
    # 使用示例
    gene_id = 'ENSG00000198062'
    transcript_dict = {
        "ENST00000343518": "GGTAGACGCGATCTGCTGGCTACTACCGGCCTTCCCTGGCTGTTAAAAGCAGATGGTGGCTGAGGCTGGTTCAATGCCGGCTGCCTCCTCTGTGAAGAAGCCATTTGGTCTCAGAAGCAAGATGGGCAAGTGGTGCCGCCACTGCTTCGCCTGGTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGATTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGGTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGATTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGCTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGACTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGCTGCAGGGGGAGCGGCAAGAACAAAGTGGGCCCTTGGGGAGACTACGACGACAGCGCTTTCATGGAGCCGAGGTACCACGTCCGTCGAGAAGATCTGGACAAGCTCCACAGAGCTGCCTGGTGGGGTAAAGTCCCCAGAAAGGATCTCATCGTCATGCTCAAGGACACTGACATGAACAAGAAGGACAAGCAAAAGAGGACTGCTCTACATCTGGCCTCTGCCAATGGAAATTCAGAAGTAGTAAAACTCCTGCTGGACAGACGATGTCAACTTAATATCCTTGACAACAAAAAGAGGACAGCTCTGACAAAGGCCGTACAATGCCAGGAAGATGAATGTGCGTTAATGTTGCTGGAACATGGCACTGATCCGAATATTCCAGATGAGTATGGAAATACCGCTCTACACTATGCTATCTACAATGAAGATAAATTAATGGCCAAAGCACTGCTCTTATACGGTGCTGATATCGAATCAAAAAACAAGCATGGCCTCACACCACTGTTACTTGGTGTACATGAGCAAAAACAGCAAGTGGTGAAATTTTTAATCAAGAAAAAAGCAAATTTAAATGCACTGGATAGATATGGAAGAACTGCTCTCATACTTGCTGTATGTTGTGGATCGGCAAGTATAGTCAGCCTTCTACTTGAGCAAAACATTGATGTATCTTCTCAAGATCTATCTGGACAGACGGCCAGAGAGTATGCTGTTTCTAGTCGTCATAATGTAATTTGCCAGTTACTTTCTGACTACAAAGAAAAACAGATACTAAAAGTCTCTTCTGAAAACAGCAATCCAGAACAAGACTTAAAGCTGACATCAGAGGAAGAGTCACAAAGGCTTAAAGGAAGTGAAAATAGCCAGCCAGAGGAAATGTCTCAAGAACCAGAAATAAATAAGGGTGGTGATAGAAAGGTTGAAGAAGAAATGAAGAAGCACGGAAGTACTCATATGGGATTCCCAGAAAACCTGACTAACGGTGCCACTGCTGACAATGGTGATGATGGATTAATTCCACCAAGGAAAAGCAGAACACCTGAAAGCCAGCAATTTCCTGACACTGAGAATGAACAGTATCACAGTGATGAACAAAATGATACTCAGAAGCAACTTTCTGAAGAACAGAACACTGGAATATTACAAGATGAGATTCTGATTCATGAAGAAAAGCAGATAGAAGTGGCTGAAAATGAATTCTGAGCTTTCTCTTAGTTATAAGAAAGAAAAAGACCTCTTGCATGAAAATAGTACGTTGCAGGAAGAAATTGTCATGCTAAGACTGGAACTAGACATAATGAAACATCAGAGCCAGCTAAGAGAAAAGAAATATTTGGAGGAAATTGAAAGTGTGGAAAAAAAGAATGATAATCTTTTAAAGGGTCTACAACTGAATGAGCTCACCATGGATGATGATACTGCCGTGCTCGTCATTGACAAC",
        "ENST00000621704": "ATGGTGGCTGAGGCTGGTTCAATGCCGGCTGCCTCCTCTGTGAAGAAGCCATTTGGTCTCAGAAGCAAGATGGGCAAGTGGTGCCGCCACTGCTTCGCCTGGTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGATTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGGTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGACTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGCTGCAGGGGGAGCGGCAAGAACAAAGTGGGCCCTTGGGGAGACTACGACGACAGCGCTTTCATGGAGCCGAGGTACCACGTCCGTCGAGAAGATCTGGACAAGCTCCACAGAGCTGCCTGGTGGGGTAAAGTCCCCAGAAAGGATCTCATCGTCATGCTCAAGGACACTGACATGAACAAGAAGGACAAGCAAAAGAGGACTGCTCTACATCTGGCCTCTGCCAATGGAAATTCAGAAGTAGTAAAACTCCTGCTGGACAGACGATGTCAACTTAATATCCTTGACAACAAAAAGAGGACAGCTCTGACAAAGGCCGTACAATGCCAGGAAGATGAATGTGCGTTAATGTTGCTGGAACATGGCACTGATCCGAATATTCCAGATGAGTATGGAAATACCGCTCTACACTATGCTATCTACAATGAAGATAAATTAATGGCCAAAGCACTGCTCTTATACGGTGCTGATATCGAATCAAAAAACAAGCATGGCCTCACACCACTGTTACTTGGTGTACATGAGCAAAAACAGCAAGTGGTGAAATTTTTAATCAAGAAAAAAGCAAATTTAAATGCACTGGATAGATATGGAAGAACTGCTCTCATACTTGCTGTATGTTGTGGATCGGCAAGTATAGTCAGCCTTCTACTTGAGCAAAACATTGATGTATCTTCTCAAGATCTATCTGGACAGACGGCCAGAGAGTATGCTGTTTCTAGTCGTCATAATGTAATTTGCCAGTTACTTTCTGACTACAAAGAAAAACAGATACTAAAAGTCTCTTCTGAAAACAGCAATCCAGAACAAGACTTAAAGCTGACATCAGAGGAAGAGTCACAAAGGCTTAAAGGAAGTGAAAATAGCCAGCCAGAGGAAATGTCTCAAGAACCAGAAATAAATAAGGGTGGTGATAGAAAGGTTGAAGAAGAAATGAAGAAGCACGGAAGTACTCATATGGGATTCCCAGAAAACCTGACTAACGGTGCCACTGCTGACAATGGTGATGATGGATTAATTCCACCAAGGAAAAGCAGAACACCTGAAAGCCAGCAATTTCCTGACACTGAGAATGAACAGTATCACAGTGATGAACAAAATGATACTCAGAAGCAACTTTCTGAAGAACAGAACACTGGAATATTACAAGATGAGATTCTGATTCATGAAGAAAAGCAGATAGAAGTGGCTGAAAATGAATTCTGAGCTTTCTCTTAGTTATAAGAAAGAAAAAGACCTCTTGCATGAAAATAGTACGTTGCAGGAAGAAATTGTCATGCTAAGACTGGAACTAGACATAATGAAACATCAGAGCCAGCTAAGAGAAAAGAAATATTTGGAGGAAATTGAAAGTGTGGAAAAAAAGAATGATAATCTTTTAAAGGGTCTACAACTGAATGAGCTCACCATGGATGATGATACTGCCGTGCTCGTCATTGACAACGGCTCTGGCATGTGCAAGGCCGGCTTTGCAGGTGACGATGCCCCCCGGGCTGTCTTCCCTTCCATCGTGGGGTGCCCCAGGCACCAGAACATGATGGGGGGCATGCGTC",
        "ENST00000452800": "CTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGGTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGATTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGCTGCAGGGGGAGCGGCAAGAGCAACGTGGGCACTTCTGGAGACCACGACGACTCTGCTATGAAGACACTCAGGAGCAAGATGGGCAAGTGGTGCTGCCACTGCTTCCCCTGCTGCAGGGGGAGCGGCAAGAACAAAGTGGGCCCTTGGGGAGACTACGACGACAGCGCTTTCATGGAGCCGAGGTACCACGTCCGTCGAGAAGATCTGGACAAGCTCCACAGAGCTGCCTGGTGGGGTAAAGTCCCCAGAAAGGATCTCATCGTCATGCTCAAGGACACTGACATGAACAAGAAGGACAAGCAAAAGAGGACTGCTCTACATCTGGCCTCTGCCAATGGAAATTCAGAAGTAGTAAAACTCCTGCTGGACAGACGATGTCAACTTAATATCCTTGACAACAAAAAGAGGACAGCTCTGACAAAGGCCGTACAATGCCAGGAAGATGAATGTGCGTTAATGTTGCTGGAACATGGCACTGATCCGAATATTCCAGATGAGTATGGAAATACCGCTCTACACTATGCTATCTACAATGAAGATAAATTAATGGCCAAAGCACTGCTCTTATACGGTGCTGATATCGAATCAAAAAACAAGattgagtctcactctgtttcccaggctggagtgcagtggtgtgatcacagctcactgcaacctcgacctcccgaagctcagatggtcctccaacctcaacctcagttgtttttttctttctttctttttatttttatttttatttttatttttttttagtagagatttggtttttgccatgtttctcaggctggtcttgaactcctgggctcaagcaattcatccccctcagcctcccaaagtgctaggattacagCATGGCCTCACACCACTGTTACTTGGTGTACATGAGCAAAAACAGCAAGTGGTGAAATTTTTAATCAAGAAAAAAGCAAATTTAAATGCACTGGATAGATATGGAAGAACTGCTCTCATACTTGCTGTATGTTGTGGATCGGCAAGTATAGTCAGCCTTCTACTTGAGCAAAACATTGATGTATCTTCTCAAGATCTATCTGGACAGACGGCCAGAGAGTATGCTGTTTCTAGTCGTCATAATGTAATTTGCCAGTTACTTTCTGACTACAAAGAAAAACAGATACTAAAAGTCTCTTCTGAAAACAGCAATCCAGAACAAGACTTAAAGCTGACATCAGAGGAAGAGTCACAAAGGCTTAAAGGAAGTGAAAATAGCCAGCCAGAGGAAATGTCTCAAGAACCAGAAATAAATAAGGGTGGTGATAGAAAGGTTGAAGAAGAAATGAAGAAGCACGGAAGTACTCATATGGGATTCCCAGAAAACCTGACTAACGGTGCCACTGCTGACAATGGTGATGATGGATTAATTCCACCAAGGAAAAGCAGAACACCTGAAAGCCAGCAATTTCCTGACACTGAGAATGAACAGTATCACAGTGATGAACAAAATGATACTCAGAAGCAACTTTCTGAAGAACAGAACACTGGAATATTACAAGATGAGATTCTGATTCATGAAGAAAAGCAGATAGAAGTGGCTGAAAATGAATTCTGAGCTTTCTCTTAGTTATAAGAAAGAAAAAGACCTCTTGCATGAAAATAGTACGTTGCAGGAAGAAATTGTCATGCTAAGACTGGAACTAGACATAATGAAACATCAGAGCCAGCTAAGAGAAAAGAAATATTTGGAGGAAATTGAAAGTGTGGAAAAAAAGAATGATAATCTTTTAAAGGGTCTACAACTGAATGAGCTCACCATGGATGATGATACTGCCGTGCTCGTCATTGACAAC"
    }

    designer = PrimerDesigner(gene_id, transcript_dict, dict(),19,
                              min_gc = 0.4, max_gc = 0.6, min_len = 18, max_len = 25,
                              min_tm = 55, max_tm = 65, min_size = 300, max_size = 500)
    designer.design_primers_for_transcript("ENST00000343518")
    designer.design_primers_for_transcript("ENST00000621704")
    designer.design_primers_for_transcript("ENST00000452800")
    # designer.visualize_primers_multi(["ENST00000343518", "ENST00000621704", "ENST00000452800"])
    # designer.visua
    # lize_primers_single("ENST00000343518")
    a = designer.get_primer_info("ENST00000343518")
    print(a)
    designer.print_primer_info("ENST00000343518")
