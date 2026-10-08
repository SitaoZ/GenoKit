Introduction
============

What GenoKit does
-----------------

GenoKit uses a reference annotation and its matching genome to support four
connected groups of operations:

* **Database and utilities:** ``create``, ``stat``, ``fasta``, ``iupac`` and ``kmer``.
* **Extraction:** ``gene``, ``transcript``, ``exon``, ``intron``, ``cds``, ``utr``,
  ``uorf``, ``dorf``, ``promoter``, ``terminator`` and ``igr``.
* **Design and sequence searches:** ``primer``, ``sgrna``, ``sirna``, ``shrna``,
  ``barcode`` and ``motif``.
* **Visualization:** ``view``, ``circos``, ``ppi`` and ``hgvs``.

There are eleven extraction commands; ``utr`` covers both untranslated ends.
``gene`` returns the genomic interval, whereas ``transcript`` reconstructs a
spliced sequence. These are different biological products, not alternative
names for the same extraction.

How the components connect
--------------------------

.. code-block:: text

   GTF/GFF3 -- create --> annotation SQLite database
   FASTA -------------> indexed sequence access
                              |
                 extraction / design / gene view
                              |
                 FASTA, CSV, feature tables, figures

The default annotation database is a GenoKit SQLite schema with lazy access to
models. It is not a gffutils FeatureDB or an R TxDb, even though those formats
also use SQLite. Large reference sequences are accessed through an indexed
FASTA reader instead of being concatenated into one in-memory genome string.
Individual commands can still retain substantial candidate/model data.

K-mer catalogues are separate databases used by optional design screening.
Motif searches operate on sequence files; barcodes are generated independently
of an annotation. Circos and PPI use supplied track/edge tables. HGVS uses
accession-specific sequence records rather than the annotation database.

Scope of predictions
--------------------

Design commands generate sequence-based **candidates**. K-mer occurrence counts
are not mismatch-aware alignment searches, and the current primer Tm estimate
is not a salt- and concentration-adjusted thermodynamic model. ORF detection
reports sequence-defined candidates, not evidence of translation. Match the
reference, annotation release and biological target before interpreting results.

A browser interface is available at `GenoKit Web <http://www.genokit.cn/>`_.
The instructions here describe the local CLI; web deployments can offer a
subset of commands and different limits.
