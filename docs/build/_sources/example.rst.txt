Worked example with a small genome
==================================

Download the inputs
-------------------

* :download:`demo.fa <_downloads/demo.fa>`: a synthetic 1,000-nt chromosome.
* :download:`demo.gtf <_downloads/demo.gtf>`: two two-exon transcripts, one on
  each strand, with CDS, stop-codon and UTR annotations.

Save both files in a new writable directory and run the commands below there.
In a source checkout, the files also reside in ``docs/source/_downloads``.
These synthetic sequences demonstrate file formats and strand handling; their
repeated composition is not intended as a realistic primer-design benchmark.

Build and extract
-----------------

.. code-block:: bash

   mkdir -p results
   GenoKit create -s gtf -a demo.gtf -g demo.fa -o results -p demo
   GenoKit transcript -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -f fasta -o results/transcripts.fa
   GenoKit gene -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -i DEMO2 -f fasta --print
   GenoKit intron -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -f fasta -o results/introns.fa
   GenoKit cds -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -f fasta -o results/cds.fa
   GenoKit uorf -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -l 6 -f csv -o results/uorf.csv
   GenoKit dorf -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -l 6 -f csv -o results/dorf.csv

The genome contains DEMO1 (101--340, positive strand) and DEMO2
(501--740, negative strand). Their transcripts DEMO1.t1 and DEMO2.t1 each have
two 90-nt exons, giving a 180-nt mature sequence; each has a 60-nt intron.
The two mature sequences were constructed to be identical after orienting the
negative-strand transcript. GTF CDS rows contain 120 nt; the separate stop codon
adds 3 nt to the extracted CDS. An exact-content check is more informative than
just checking that an output file exists.

Explore the remaining features
------------------------------

.. code-block:: bash

   GenoKit utr -d results/demo.gtf.sqlite -g demo.fa -s gtf -f csv -o results/utr.csv
   GenoKit promoter -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -l 50 -u 0 -f fasta -o results/promoters.fa
   GenoKit terminator -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -l 50 -u 0 -f fasta -o results/terminators.fa
   GenoKit igr -d results/demo.gtf.sqlite -g demo.fa -s gtf -f fasta -o results/igr.fa
   GenoKit motif -i results/transcripts.fa -m ATG -t dna -o results/motif.csv
   GenoKit view -d results/demo.gtf.sqlite -g demo.fa -s gtf \
     -i DEMO1 -f pdf -o results/view

The intergenic gap is 341--500 (160 nt). Chromosome-terminal sequence is not
included by the IGR command. Upstream/downstream windows are relative to gene
strand, not always lower/higher chromosome coordinates.

A small k-mer catalogue
-----------------------

.. code-block:: bash

   GenoKit kmer -f demo.fa --source-type genome -l 20 -r \
     --engine python -o results/kmer
   GenoKit kmer -f results/transcripts.fa --source-type transcriptome -l 20 -r \
     --engine python -o results/kmer

This registers two biological sources at k=20 in one catalogue. It does not
provide k=18, 19 or 21--25, so a primer run requesting those lengths would need
additional datasets or restrict candidate length to 20.
