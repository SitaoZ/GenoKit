Primer design
=============

Select the actual template
--------------------------

.. list-table::
   :header-rows: 1

   * - ``--template-type``
     - Template
     - ``--kmer-source``
   * - genomic
     - Continuous, strand-oriented gene interval including introns
     - genome
   * - cdna (default)
     - Each mature, exon-joined transcript of the gene
     - transcriptome

If omitted, the k-mer source is inferred from the template type. An explicitly
conflicting source is rejected. Choosing cDNA does not require that a primer
cross an exon junction, and does not guarantee absence of genomic amplification.

.. code-block:: bash

   GenoKit primer -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID --template-type cdna --kmer-source transcriptome \
     -k results/kmer/genokit.kmer.sqlite -n 18 -m 25 -x 100 -y 200 \
     -f csv -o results/primer.csv
   GenoKit primer -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID --template-type genomic --kmer-source genome \
     -k results/kmer/genokit.kmer.sqlite -f csv -o results/genomic_primer.csv

``-k`` is optional. Without it, design still operates but catalogue frequencies
and off-target counts are unavailable. Use ``--print`` instead of ``-o`` for
stdout, never both. Selecting a gene also generates a primer plot in the working
directory; the plot filename is based on the gene ID.

Candidate search and ranking
----------------------------

The CLI searches the full oriented template with sliding windows. It is not
limited to transcript ends, and an internal amplicon can be found in a long
transcript. The Python method's optional ``target_region`` is not a CLI flag.

The default CLI constraints are 18--25 nt primer length, GC fraction 0.4--0.6,
Tm 50--70 degrees C, and product size 100--1000 bp. Tm uses the Wallace rule:

.. math::

   T_m = 2(A+T) + 4(G+C).

Non-ATCG windows and runs of three identical bases are rejected. A simple
repeated-4-mer heuristic called ``fast_hairpin`` in the source rejects additional
candidates; it is not a thermodynamic hairpin or primer-dimer calculation.
When a catalogue is supplied, candidate lengths must have exact-k support and
candidates must have a positive occurrence count.

Pairs must have nonoverlapping, inward-facing sites and a product length within
the requested range. The score is

.. math::

   S = 2|T_{m,F}-T_{m,R}| + 0.1|GC_F-GC_R|,

where GC is the rounded fraction used by the implementation. The lowest score
wins. **One pair per template** is retained; ties do not replace the earlier
pair. K-mer counts and coverage annotations are not terms in this score. A plot
can therefore show an end-biased chosen pair even though internal candidates
were searched. There is no CLI option to output every tied pair.

Coordinates and core output fields
----------------------------------

CSV/FASTA ``Start`` and ``End`` are **1-based inclusive template coordinates**.
The conversion from internal offsets occurs only at serialization. For both
primer directions, ``Length = End - Start + 1``; sequences are written 5-prime
to 3-prime. ``Fragment`` is the complete amplicon length including both primers
and, in genomic mode, intervening introns. ``Template Type`` identifies the
mode; genomic rows have no transcript ID. ``Template Seq`` is the full template,
so CSV files can be large when many templates are reported.

Interpret exact-match specificity fields
----------------------------------------

* ``Kmer Freq``: total canonical occurrences in the selected reference.
* ``Target Template Count``: the gene's transcript templates in cDNA mode,
  or its single genomic interval in genomic mode.
* ``Target Templates Hit`` and ``Target Coverage``: templates containing at
  least one exact occurrence in either orientation, and their fraction.
* ``Target Exact Matches``: all exact occurrences within the target templates,
  allowing more than one match per template.
* ``Off-target Exact Matches``: catalogue frequency minus target exact matches,
  only when the source-aware canonical counting contract is available.
* ``Target Pair Templates Hit`` and ``Target Pair Coverage``: templates in which
  the same pair has inward-facing, nonoverlapping matches satisfying product
  size limits; coverage is a fraction, not a percentage.

``Specificity Status`` is ``ok`` for compatible subtraction, or explains an
unavailable result: ``no_database``, ``unknown_counting``,
``incompatible_counting``, ``missing_frequency`` or ``count_mismatch``.
``NA`` is not zero. A count mismatch can indicate inconsistent references.
Nonnegative subtraction alone does not prove that references are identical.

An ideal transcriptome frequency need not equal the number of isoforms: a
primer can be present only in a subset or occur repeatedly in one transcript.
Use coverage and exact-match fields together. Pair coverage describes targets;
it does not enumerate genome-wide off-target PCR products. Near matches and
mismatch-tolerant amplification are outside this counting model.

Troubleshooting no candidates
-----------------------------

Check exact-k availability, source type, canonical counting, template length,
GC/Tm constraints and product limits. Repeated or ambiguous sequence can be
filtered out. Reducing the minimum product length cannot fix a missing k value.
``--specificity-mode`` and ``--seed-size`` are retired options; no extra seed
index is required by this workflow.
