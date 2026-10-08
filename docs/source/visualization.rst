Gene, genome and network visualization
======================================

Gene structures
---------------

.. code-block:: bash

   GenoKit view -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i GENE_ID -l 6 --orf-label-limit 20 -f pdf -o results/gene_view

``-o`` is a directory. Formats are PDF (default), PNG and TIFF. The plot shows
transcript/exon structures, coding and untranslated regions, and candidate
regulatory ORFs. ``--orf-label-limit 0`` suppresses individual ORF numbers,
not the underlying features. Omitting the gene selector processes all genes
and can create many files; begin with a selected gene when checking inputs.

Circular genome plots
---------------------

.. code-block:: bash

   GenoKit circos -a example/chromosome_length.csv \
     -b example/barplot_data.csv -m example/heatmap_data.csv \
     -k example/link_data.csv -f pdf -o results/circos.pdf

The ``example/...`` paths first resolve in the working directory, then in the
package's example directory. Inputs are comma-separated tables with a header.

.. list-table::
   :header-rows: 1

   * - Option
     - Columns, in order
   * - ``-a/--chrom``
     - chr,start,end (end supplies chromosome length)
   * - ``-c/--cytoband``
     - chromosome,start,end,band name,stain
   * - ``-b/--bar``, ``-p/--point``, ``-l/--line``, ``-m/--heatmap``
     - chr,start,end,value
   * - ``-k/--link``
     - chr1,start1,end1,chr2,start2,end2

``-t/--track`` detects a line or link input from its header. Chromosome names in
tracks must match the chromosome table. Supply a filename whose extension
matches the desired format. This is a plotting workflow, not automatic
whole-genome annotation summarization from a SQLite database.

Protein interaction networks
----------------------------

.. code-block:: text

   TP53,MDM2
   TP53,ATM
   MDM2,RPL11

Save these pairs in ``interactions.csv`` **without a header**:

.. code-block:: bash

   GenoKit ppi -i interactions.csv -s comma -l auto -o results/ppi

Use ``-s tab`` for tab-separated edges. GenoKit uses supplied interactions;
it does not automatically retrieve a STRING network for these IDs. Record the
source and any filtering of the edge list separately.

Layouts include auto, fr, kk, forceatlas2, circle, star, tree, grid, shell,
spectral, spiral and random. Outputs include network and centrality PDFs,
a centrality table, a GEXF network and a log. ``ppi_centrality.csv`` currently
contains tab-separated values despite its extension. Layout coordinates can
vary with algorithms and random initialization; retain exported results.

On a headless server, use ``export MPLBACKEND=Agg`` and a writable
``MPLCONFIGDIR``. This does not affect the sequence coordinate conventions.
