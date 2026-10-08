Genomic feature extraction
==========================

Choose a biological product
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 18 64

   * - Command
     - ``-i`` selects
     - Product
   * - gene
     - Gene ID
     - Continuous genomic gene interval, strand oriented.
   * - transcript
     - Transcript ID
     - Exons joined to reconstruct a mature transcript.
   * - exon
     - Transcript ID
     - Individual exons associated with a transcript.
   * - intron
     - Transcript ID
     - Intervals between consecutive transcript exons.
   * - cds
     - Transcript ID
     - Joined annotated coding segments; see stop-codon notes below.
   * - utr
     - Transcript ID
     - Both 5-prime and 3-prime untranslated sequences.
   * - promoter
     - Gene ID
     - Window around the strand-aware gene 5-prime boundary.
   * - terminator
     - Gene ID
     - Window around the strand-aware gene 3-prime boundary.
   * - igr
     - Not available
     - Gaps between neighboring genes on the same chromosome.
   * - uorf / dorf
     - Transcript ID
     - Candidate ORFs relative to the canonical CDS; see :doc:`regulatory_orfs`.

Omitting ``-i`` processes applicable records across the database. Shared exons
can occur in several transcript-associated records; this is not a promise of
genome-wide deduplication. Use the same inclusion and deduplication rules when
comparing counts with other programs.

Core examples
-------------

.. code-block:: bash

   GenoKit gene -d results/human.gtf.sqlite -g human.fa -s gtf -f csv -o results/gene.csv
   GenoKit transcript -d results/human.gtf.sqlite -g human.fa -s gtf -f fasta -o results/transcript.fa
   GenoKit exon -d results/human.gtf.sqlite -g human.fa -s gtf -f fasta -o results/exon.fa
   GenoKit intron -d results/human.gtf.sqlite -g human.fa -s gtf -f fasta -o results/intron.fa
   GenoKit cds -d results/human.gtf.sqlite -g human.fa -s gtf -f fasta -o results/cds.fa
   GenoKit utr -d results/human.gtf.sqlite -g human.fa -s gtf -f csv -o results/utr.csv

``transcript --upper`` (``-u``) marks CDS bases in uppercase and UTR bases in
lowercase. It does not remove UTRs. Use ``cds`` to request coding sequence alone.
The UTR command returns both sides; it has no ``--utr5-only`` switch.

CDS reconstruction
------------------

CDS extraction concatenates annotated segments and handles separate GTF
``stop_codon`` features. The resulting sequence can therefore differ from a
comparison program that concatenates only CDS rows. Annotation phase is retained
in feature records; do not assume that the extraction routine automatically
repairs a partial or inconsistent reading frame. RNA-design CDS templates do
not separately append GTF stop-codon rows; see :doc:`rna_design`.

Promoter and terminator boundaries
----------------------------------

.. code-block:: bash

   GenoKit promoter -d results/human.gtf.sqlite -g human.fa -s gtf \
     -l 1000 -u 0 -f fasta -o results/promoter.fa
   GenoKit terminator -d results/human.gtf.sqlite -g human.fa -s gtf \
     -l 1000 -u 0 -f fasta -o results/terminator.fa

Both commands default to 100 nt outside and 10 nt inside the gene boundary.
For promoter, ``-l`` extends upstream and ``-u`` includes the beginning of the
gene; for terminator, ``-l`` extends downstream and ``-u`` includes its end.
``-u 0`` requests only flanking sequence. These are gene-level windows, not
one promoter per alternative transcript start. Intervals are clipped to
reference boundaries, so output can be shorter than the requested length.

Intergenic regions
------------------

.. code-block:: bash

   GenoKit igr -d results/human.gtf.sqlite -g human.fa -s gtf -f fasta -o results/igr.fa

The current implementation sorts genes by start and considers consecutive
pairs, omitting overlapping/adjacent pairs and chromosome-terminal gaps.
It does not first merge all overlapping gene intervals. Nested genes can
therefore make its output differ from a complement of the union of all genes.
IGR sequences follow reference orientation; strand metadata identifies the
two flanking gene strands and can contain a value such as ``+|-``. Such a
paired value is not a standard single GFF strand character.

.. important::

   The parser advertises ``igr -l/--igr_length``, but the current ``get_igr``
   extraction path does not apply it. Do not rely on that flag to filter lengths;
   filter the resulting sequences explicitly if a threshold is required.
