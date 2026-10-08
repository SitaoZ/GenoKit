HGVS variant-context visualization
==================================

Inputs and reference records
----------------------------

Use an accession-qualified HGVS expression with ``-s/--string`` or a file with
``-i``. A file can contain one expression per line; CSV/TSV readers also accept
the first column and recognized hgvs/variant/variation headers.

.. code-block:: bash

   GenoKit hgvs -i variants.txt -o results/hgvs
   GenoKit hgvs -i variants.txt --genbank reference.gb -o results/hgvs_local

The reference must match the accession, its version and the coordinate type.
A local GenBank/GenPept record avoids a reference download. Otherwise records
are retrieved from NCBI and cached in ``.genokit_cache/hgvs``; use ``--cache-dir``
to choose another location. Configure requests with ``--email``/``NCBI_EMAIL``
and, if available, ``--api-key``/``NCBI_API_KEY``. Do not publish API keys in
example scripts or archived logs.

Supported contexts
------------------

The parser recognizes DNA coordinate types g., c., n., m., o., RNA r. and
protein p. The available forms depend on coordinate type; this does not imply
complete HGVS grammar coverage. ``--variation`` can validate an expected kind:
sub, del, ins, delins, dup, inv, rpt, alleles, ext or fs.

Output supports PNG, PDF and SVG. For multiple inputs, use an output directory
or prefix. ``--context`` defaults to 20 and ``--dpi`` to 200. Local sequence,
codon and amino-acid context depend on the reference record and variant type.

Reference validation
--------------------

A mismatch can indicate the wrong accession version, coordinate type or
reference allele. Fix those inputs first. ``--allow-reference-mismatch`` and
``--skip-three-prime-check`` are diagnostic overrides; record their use and
do not interpret an overridden check as successful reference validation.
This visualization command is separate from the legacy GenoKitGB entry point.
