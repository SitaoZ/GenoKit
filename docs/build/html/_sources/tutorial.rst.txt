Tutorial: annotation to sequence
================================

Prepare matching inputs
-----------------------

Use a GTF/GFF3 annotation and reference FASTA from the same assembly and release.
The chromosome identifiers must match exactly (``1`` and ``chr1`` are different
keys). Create a writable output directory. The paths in this tutorial are
placeholders; :doc:`example` provides actual small downloadable inputs.

.. code-block:: bash

   mkdir -p results
   GenoKit create -s gtf -a human.gtf -g human.fa -o results -p human
   GenoKit stat -d results/human.gtf.sqlite -g human.fa -s gtf -o results/stat.csv

The first command creates ``results/human.gtf.sqlite`` and checks/builds the
FASTA index. Reuse this database when the annotation is unchanged. Rebuild it
when changing the annotation; do not simply rename a database from another tool.

Extract sequences
-----------------

.. code-block:: bash

   GenoKit transcript -d results/human.gtf.sqlite -g human.fa -s gtf \
     -f fasta -o results/transcripts.fa
   GenoKit cds -d results/human.gtf.sqlite -g human.fa -s gtf \
     -f fasta -o results/cds.fa
   GenoKit promoter -d results/human.gtf.sqlite -g human.fa -s gtf \
     -l 1000 -u 0 -f fasta -o results/promoters.fa

The transcript output includes spliced exons; the CDS output uses coding
segments and the implemented stop-codon handling. The promoter example requests
upstream sequence only. Parent directories are not generally created by every
output writer, so create them explicitly.

Select one record
-----------------

.. code-block:: bash

   GenoKit gene -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID -f fasta --print
   GenoKit transcript -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i TRANSCRIPT_ID -f fasta --print

Replace these identifiers with exact IDs from the annotation, including version
suffixes and any prefixes. ``-i`` is a gene selector for ``gene``, promoter,
terminator and gene-based design commands, but a transcript selector for
transcript, exon, intron, CDS, UTR and ORF extraction.

Continue to analysis and design
-------------------------------

Search the extracted transcript FASTA directly:

.. code-block:: bash

   GenoKit motif -i results/transcripts.fa -m ATG -t dna -o results/atg_sites.csv

For primer design, choose :doc:`primer`'s genomic or cDNA template. For optional
occurrence screening, first build the matching :doc:`kmer` source and k values.
To draw a selected model, use :doc:`visualization`. Do not infer that an empty
design output means the gene was missing: valid models can have no candidates
passing all filters.
