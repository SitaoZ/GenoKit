Troubleshooting and frequently asked questions
==============================================

Why are old commands rejected?
------------------------------

Commands and annotation styles are lowercase. Use ``cds``, ``utr``, ``uorf``,
``dorf`` and ``igr``, not historical mixed-case spellings. ``transcript`` replaces
old cDNA/mRNA examples. For extraction, ``-g`` is reference FASTA and ``-f`` is
output format. See :doc:`migration` and the installed command's help.

Why does an ID not exist?
-------------------------

Check whether the command requires a gene or transcript ID, and whether the
annotation includes prefixes or version suffixes. A symbol such as TP53 is not
automatically resolved to its annotation ID. Confirm that the database was built
from the intended annotation rather than another release.

Why is FASTA indexing failing?
------------------------------

Check contig names, file access, writable index paths and compression. Ordinary
gzip FASTA is not BGZF. Use plain FASTA or bgzip, and rebuild stale indexes after
changing the reference. A chromosome alias is not automatically reconciled.

Why is an SQLite database rejected?
-----------------------------------

GenoKit annotation SQLite, GenoKit k-mer SQLite, gffutils SQLite and R TxDb
SQLite have different schemas. Use ``create`` for annotation databases and
``kmer`` for count catalogues. File renaming does not convert schemas.

Why are no primers or RNA targets produced?
-------------------------------------------

Check target length, region, template type and all sequence filters. Primer
catalogues must support exact candidate lengths, the matching source and
canonical counts. For siRNA/shRNA, all-mode may have no shared candidate across
isoforms; any-mode changes the biological coverage requirement. An empty result
can be legitimate, rather than a reason to weaken filters without inspection.

Why do primers appear near transcript ends?
-------------------------------------------

The CLI searches the full template, but retains one lowest-scoring pair per
template, with earlier ties retained. This does not imply that only endpoints
were searched. See :doc:`primer` for the score and coverage fields.

Why are off-target counts NA?
-----------------------------

Read ``Specificity Status``. No database, unknown canonical/source metadata,
missing frequency or counts inconsistent with target occurrences prevent valid
subtraction. NA is not evidence of zero off-target occurrences. Raw legacy count
files do not provide all the metadata of a source-aware catalogue.

Why are outputs shorter or different from another tool?
-------------------------------------------------------

Inspect chromosome-edge clipping, transcript versus genomic coordinates,
CDS stop-codon rules, UTR annotations, isoform-specific exons and IGR overlap
handling. Some CSV coordinates describe the parent transcript span, not every
CDS block. Use GFF/GTF block coordinates where appropriate.

Why is a figure/log in the working directory?
---------------------------------------------

Several commands create logs or automatically named plots independently of
``-o``. Run separate jobs in separate directories. For headless rendering:

.. code-block:: bash

   mkdir -p .mplconfig
   export MPLBACKEND=Agg
   export MPLCONFIGDIR="$PWD/.mplconfig"

Can I use GenoKitGB?
--------------------

The distribution declares a ``GenoKitGB`` entry point, but the current legacy
implementation imports obsolete module paths. It is not documented here as a
working feature-extraction workflow. ``GenoKit hgvs --genbank`` is a separate,
active reference-input option.

What should an issue report include?
------------------------------------

Include the exact command, full traceback, GenoKit/Python versions, OS, relevant
input format and a minimal non-sensitive example. Report whether an empty output
or an error occurred, and retain the command's log if one was created.
