Motif searches and barcodes
===========================

Search known patterns
---------------------

.. code-block:: bash

   GenoKit motif -i sequences.fa -m ATG -t dna -o results/atg.csv
   GenoKit motif -i sequences.fa -m ATNG -d -b -t dna -o results/degenerate.csv
   GenoKit motif -i sequences.fa -f patterns.txt -x 1 -o results/patterns.csv
   GenoKit motif -i sequences.fa -m 'ATG[ACGT]{3}TAA' -e -o results/regex.csv

Use ``-m`` for one pattern or ``-f`` for a text file containing one pattern per
line. ``-b`` searches both strands, ``-c`` requests case-sensitive matching,
``-d`` enables IUPAC ambiguity codes, ``-e`` enables regular expressions, and
``-x`` sets the allowed mismatch count (default 0). ``-t rna`` selects RNA mode.
Use a single matching mode deliberately; do not assume regular expressions and
mismatch matching have identical semantics. Quote expressions to prevent shell
expansion. This is supplied-pattern search, not de novo motif discovery or a
PWM database-scanning interface.

Generate barcode candidates
---------------------------

.. code-block:: bash

   GenoKit barcode -l 8 -n 20 -s dna -f csv -o results/barcodes.csv
   GenoKit barcode -l 10 -n 20 -s rna -f fasta -o results/rna_barcodes.fa

Supply ``-l`` and ``-s dna|rna`` explicitly; count defaults to 20. The CLI uses
GC bounds 0.4--0.6 and maximum homopolymer length 2. Randomly generated accepted
barcodes are unique within the generated list. The generator limits attempts,
so demanding requests can yield fewer records than requested.

The CLI does not expose all methods of the underlying BarcodeDesigner class,
such as every paired-index or distance-selection helper. Do not assume that a
random barcode set guarantees a specified pairwise Hamming distance. The CLI
has no seed flag; retain the actual output when reproducible barcode identities
are required.
