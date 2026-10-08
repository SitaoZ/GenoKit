Migration from older examples
=============================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Older instruction
     - Current interface
   * - ``UTR``, ``CDS``, ``uORF``, ``dORF``, ``IGR``
     - Lowercase ``utr``, ``cds``, ``uorf``, ``dorf``, ``igr``.
   * - ``cdna`` or ``mrna`` extraction
     - ``transcript``.
   * - ``vision``
     - ``view``.
   * - ``create -f GTF -g annotation.gtf``
     - ``create -s gtf -a annotation.gtf -g genome.fa -o results -p sample``.
   * - Extraction ``-f genome.fa`` / ``-t fasta``
     - ``-g genome.fa -f fasta``.
   * - gffutils-built annotation database
     - GenoKit ``create`` builds its own SQLite schema.
   * - New database named ``.bin``
     - Default ``PREFIX.STYLE.sqlite``; optional pickle ``.pkl``.
   * - ``primer --specificity-mode`` / ``--seed-size``
     - ``--template-type genomic|cdna`` and matching ``--kmer-source``;
       exact-length counts instead of the retired seed workflow.
   * - Primer output Start/End assumed 0-based
     - Current serialized coordinates are 1-based inclusive on the template.
   * - ``shrna --all`` / ``--common``
     - ``--target-mode all|any``; any is not a majority filter.
   * - shRNA full-transcript default
     - Explicit ``--region transcript``; current default is cds.
   * - Boolean uORF/dORF ``-m`` / ``-n``
     - Supply a prefix string after either option.
   * - ``--process`` / ``--rna_feature`` on extraction
     - Not active extraction options in this source release.

Short flags are command-specific: ``-p`` means output prefix for create,
stdout for gene, PAM pattern for sgrna and loop sequence for shrna.
Prefer descriptive long options when porting scripts and verify with
``GenoKit COMMAND -h``. Do not change a biological coverage rule solely to
make an old command parse.
