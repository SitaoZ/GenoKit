sgRNA, siRNA and shRNA design
=============================

sgRNA candidates
----------------

.. code-block:: bash

   GenoKit sgrna --pam
   GenoKit sgrna -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID -p NGG -l 20 -k results/kmer/genokit.kmer.sqlite \
     -f csv -o results/sgrna.csv

``--pam`` displays a PAM reference table. ``-p`` selects the PAM pattern
(default NGG); ``-l`` selects guide length (default 20). Candidate scanning
considers strand orientation. Use genome k-mer counts of the guide length for
optional occurrence screening. This is not a genome-wide mismatch-aware
alignment search or experimental efficiency validation.

Select an RNA target region
---------------------------

Both siRNA and shRNA support ``--region cds|transcript`` (default ``cds``).
CDS mode joins annotated CDS segments; a transcript without CDS falls back to
its exons. Transcript mode joins all exons, including UTRs. These RNA-design
CDS sequences do not separately append GTF stop-codon records.

``--target-mode all`` (default) requires each retained target sequence to occur
in every input transcript region of the gene. ``any`` requires at least one,
not a majority. Output lists actual hit transcripts. Constraints apply per gene
both with and without ``-i``. For siRNA, a short/unscannable transcript remains
in the all-mode coverage set; ambiguous windows are skipped rather than bases
being deleted and artificial junctions introduced.

siRNA
-----

.. code-block:: bash

   GenoKit sirna -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID -l 21 --region transcript --target-mode all --top 10 \
     -k results/kmer/genokit.kmer.sqlite -f csv -o results/sirna.csv

Allowed target lengths are 19, 20 and 21 nt, default 21. ``--top`` limits reported
candidates per gene. Optional k-mer counts should come from the transcriptome
and match the target length. Candidate scores are implementation-defined
sequence heuristics, not measured knockdown efficiencies.

shRNA
-----

.. code-block:: bash

   GenoKit shrna -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID -l 21 --region cds --target-mode any -p CTCGAG \
     -k results/kmer/genokit.kmer.sqlite -f csv -o results/shrna.csv

The length choices and target-mode interface match siRNA. ``-p/--loop`` is a
loop/linker sequence (default CTCGAG), not a PAM or thread count. Named loop
presets include ``standard``, ``miR30`` and ``simple``; custom ATCG strings are
accepted. Output includes the selected target, loop, and sense/antisense
oligonucleotide constructs. Inspect the reported construct before ordering.

All three commands accept ``--jellyfish PATH`` when querying .jf datasets.
K-mer databases are optional. To compare all versus any, keep the annotation,
region, candidate length and other settings identical.

For old shRNA workflows, ``--region transcript`` retains full-transcript target
selection. Replace retired ``--all``/``--common`` flags with the explicit
``--target-mode`` interface; ``any`` is broader than the old majority rule.
