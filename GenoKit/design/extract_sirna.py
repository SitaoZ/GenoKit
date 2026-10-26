# -*- coding: utf-8 -*-
"""Command implementation for ``GenoKit sirna``."""

import csv
import logging
import sys
from collections import defaultdict
from typing import Dict, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm

from GenoKit.database.database import genome_dict
from GenoKit.database.annotation_db import get_gene_symbol, load_annotation_database
from GenoKit.database.kmer_database import open_kmer_database
from GenoKit.design.design_sirna import SirnaDesigner


_CSV_HEADER = [
    "Gene ID", "Gene Symbol", "Transcript ID", "Target Sequence", "Guide Strand", "Passenger Strand",
    "GC", "Target Position", "Targets", "Priority", "Kmer Freq",
]


def parse_transcript(params: Tuple) -> Dict[str, Dict]:
    """Build mature, 5'-to-3' transcript sequences from exon annotations."""
    gene_isoform = defaultdict(dict)
    gene_id, transcript_ids, database, genome, _style, _output_format = params
    for transcript_id in transcript_ids:
        mrna = database[gene_id][transcript_id]["mRNA"]
        sequence = "".join(
            str(genome[exon.chr][exon.start - 1 : exon.end])
            for exon in database[gene_id][transcript_id]["exon"]
        )
        if mrna.strand == "-":
            sequence = str(Seq(sequence).reverse_complement())
        gene_isoform[gene_id][transcript_id] = sequence.upper()
    return gene_isoform


def parse_cds(params: Tuple) -> Dict[str, Dict]:
    """Build CDS sequences; use exons for non-coding transcripts as before."""
    gene_isoform = defaultdict(dict)
    gene_id, transcript_ids, database, genome, _style, _output_format = params
    for transcript_id in transcript_ids:
        mrna = database[gene_id][transcript_id]["mRNA"]
        features = database[gene_id][transcript_id].get("CDS") or database[gene_id][transcript_id]["exon"]
        sequence = "".join(
            str(genome[feature.chr][feature.start - 1 : feature.end])
            for feature in features
        )
        if mrna.strand == "-":
            sequence = str(Seq(sequence).reverse_complement())
        gene_isoform[gene_id][transcript_id] = sequence.upper()
    return gene_isoform


def _format_positions(positions: Dict[str, list[int]]) -> str:
    return "|".join(
        f"{transcript_id}:{','.join(str(position) for position in starts)}"
        for transcript_id, starts in positions.items()
    )


def _records_for_gene(
    gene_id, gene_symbol, gene, transcript_ids, transcripts, args, kmer_counts,
):
    designer = SirnaDesigner(
        gene_id, transcripts, kmer_bin=kmer_counts, si_length=args.length,
    )
    designs = designer.design(target_mode=args.target_mode)
    if args.top is not None:
        designs = designs[: args.top]

    for index, design in enumerate(designs, start=1):
        target = design["sequence"]
        guide = str(Seq(target).reverse_complement())
        if args.output_format == "fasta":
            description = (
                f"target:{target} guide:{guide} gene_strand:{gene.strand} "
                f"transcript:{'|'.join(sorted(design['positions']))} "
                f"positions:{_format_positions(design['positions'])} "
                f"gc:{design['gc_content']:.4f} priority:{design['priority']} "
                f"kmer_freq:{design['kmer_freq']}"
            )
            yield SeqRecord(Seq(guide), id=f"{gene_id}_{index}", description=description)
        else:
            yield [
                gene_id,
                gene_symbol,
                "|".join(sorted(design["positions"])),
                target,
                guide,
                target,
                design["gc_content"],
                _format_positions(design["positions"]),
                design["targets"],
                design["priority"],
                design["kmer_freq"],
            ]


def _write_records(args, records, csv_writer=None, fasta_handle=None):
    if args.output_format == "csv":
        csv_writer.writerows(records)
    else:
        SeqIO.write(records, fasta_handle, "fasta")


def get_sirna(args):
    """Design siRNAs and write records as they are generated."""
    # logging
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        filename="siRNA.log"
    )

    # logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    database, _t2g = load_annotation_database(args.database)
    if args.gene and args.gene not in database:
        raise ValueError(f"Gene {args.gene!r} is not in the annotation database.")

    genome = genome_dict(args.genome)
    # siRNAs bind mature mRNA, not unspliced genomic DNA.
    kmer_counts = open_kmer_database(
        args.kmer_bin,
        k=args.length,
        expected_source="transcriptome",
        jellyfish=getattr(args, "jellyfish", None),
    )
    sequence_parser = parse_transcript if args.region == "transcript" else parse_cds
    gene_ids = [args.gene] if args.gene else database.keys()

    if args.print:
        output_handle = sys.stdout
        close_handle = False
    elif args.output:
        output_handle = open(args.output, "w", newline="", encoding="utf-8")
        close_handle = True
    else:
        raise ValueError("-o or -v should be specified")

    try:
        csv_writer = csv.writer(output_handle) if args.output_format == "csv" else None
        if csv_writer:
            csv_writer.writerow(_CSV_HEADER)
        fasta_handle = output_handle if args.output_format == "fasta" else None

        iterator = gene_ids
        if not args.gene:
            iterator = tqdm(gene_ids, ncols=80, total=len(database), desc="siRNA processing:", colour="GREEN")
        for gene_id in iterator:
            gene = database[gene_id]["gene"]
            symbol = get_gene_symbol(database, gene_id)
            transcript_ids = [
                key for key in database[gene_id]
                if key not in {"gene", "symbol"}
            ]
            params = (gene_id, transcript_ids, database, genome, args.style, args.output_format)
            transcripts = sequence_parser(params)[gene_id]
            if transcripts:
                _write_records(
                    args,
                    _records_for_gene(
                        gene_id, symbol, gene, transcript_ids, transcripts,
                        args, kmer_counts,
                    ),
                    csv_writer=csv_writer,
                    fasta_handle=fasta_handle,
                )
    finally:
        if close_handle:
            output_handle.close()
