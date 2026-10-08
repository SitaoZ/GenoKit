Inputs, outputs and coordinates
===============================

Annotation and FASTA
--------------------

GTF attributes identify ``gene_id`` and ``transcript_id``. GFF3 uses feature
``ID`` and ``Parent`` relationships. Supply ``-s gtf`` or ``-s gff`` explicitly
where offered. Files ending in ``.gff3`` still use style ``gff``. Complete parent
relationships and exon/CDS structures are important for transcript-aware work.
A file extension alone is not a guarantee that its content is valid.

Annotation input to ``create`` can be plain or gzip-compressed. Reference FASTA
can be plain or BGZF-compressed. Ordinary gzip FASTA does not support the same
indexed access. Keep the ``.fai`` index and, for BGZF, the ``.gzi`` file beside
the reference. Regenerate indexes if reference bytes or wrapping change.

Output selection
----------------

Extraction parsers accept ``csv``, ``fasta``, ``gff`` and ``gtf``. Choose ``-f``
explicitly; the gene and terminator parsers do not provide a format default.
CSV exposes feature-specific metadata and sequences, FASTA is suitable for
sequence tools, and GFF/GTF represents genomic feature intervals rather than
an embedded FASTA sequence. Read column names rather than assuming a universal
schema across commands.

Most design commands provide CSV/FASTA. ``view`` and ``circos`` advertise PDF,
PNG and TIFF; HGVS also supports SVG. Some commands write logs or gene-named
plots in the current working directory. Use separate working directories for
independent runs and provide a new output filename to avoid overwriting results.

Coordinate systems
------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Context
     - Interpretation
   * - Input GTF/GFF3; genomic feature intervals
     - 1-based, inclusive: length = end - start + 1.
   * - Indexed Python sequence slices
     - Internal 0-based, half-open slices: ``genome[start - 1:end]``.
   * - Primer CSV/FASTA Start and End
     - 1-based, inclusive positions on the oriented design template;
       not chromosome coordinates. Internal primer offsets remain 0-based.
   * - uORF/dORF Start and End transcript fields
     - 1-based, inclusive positions in the spliced transcript.
   * - HGVS notation
     - Follows the supplied accession and coordinate type; for example,
       coding coordinates do not start at the transcript's first nucleotide.

For a primer with Start=21 and End=40, the length is 20 nt. A reverse primer's
coordinates also ascend along the oriented template; its reported sequence is
the reverse complement of that interval. For negative-strand genes, increasing
template coordinates run opposite to increasing genomic coordinates.

Do not apply one conversion rule to every descriptive field. Some legacy FASTA
metadata fields (including CDS interval summaries) are calculated differently
from genomic output columns. dORF's ``CDS Interval`` string uses legacy 0-based
boundary values, although its dORF start/end transcript fields are 1-based.
Use feature GFF/GTF intervals when genomic segment coordinates are required;
do not interpret transcript spans as continuous genomic coding intervals.

Splicing and strand
-------------------

Exons/CDS segments are read from the reference and joined in transcript
orientation, including reverse complementation on the negative strand.
A mature transcript sequence does not contain introns. A genomic gene template
does. A single spliced interval can map to multiple genomic blocks, so it cannot
in general be converted by adding only the gene start.

Identifiers and empty results
-----------------------------

Gene symbols are descriptive metadata, not guaranteed unique lookup keys.
No CDS/UTR/ORF output is expected for some noncoding or incomplete models.
Single-exon transcripts have no introns. Exact feature inclusion depends on the
command and annotation; see :doc:`extraction` and :doc:`regulatory_orfs`.
