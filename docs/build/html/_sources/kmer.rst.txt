K-mer catalogues and portability
================================

Separate the biological sources
-------------------------------

``genome`` counts reference genomic sequence. ``transcriptome`` counts spliced
transcript sequences and can count the same shared sequence in several isoforms.
Use genome counts for genomic primers and sgRNAs; use transcriptome counts for
cDNA primers, siRNAs and shRNAs. The annotation SQLite and k-mer SQLite are
separate files even when both have a ``.sqlite`` extension.

Build exact lengths required by design
--------------------------------------

.. code-block:: bash

   GenoKit kmer -f human.fa --source-type genome --lengths 18-25 \
     -o results/kmer -r --engine jellyfish -j jellyfish -t 4
   GenoKit kmer --source-type transcriptome -d results/human.gtf.sqlite -g human.fa \
     --lengths 18-25 -o results/kmer -r --engine jellyfish -j jellyfish -t 4

One ``results/kmer/genokit.kmer.sqlite`` catalogue can register multiple sources
and k values. ``--lengths`` accepts ranges and comma-separated values, whereas
``-l`` specifies one length (default 18). ``-r`` enables canonical reverse-
complement counting, required by the primer catalogue workflow.

For a prebuilt mature-transcript FASTA, use ``-f transcripts.fa`` with
``--source-type transcriptome`` instead of ``-d``/``-g``. Using genomic sequence
with a transcriptome label does not make it a transcriptome reference.

Engines
-------

``--engine auto`` selects an available implementation; explicit choices are
``jellyfish`` and ``python``. The Python implementation is limited to small
inputs up to 250 MiB and k values 1--31. It is useful for demonstrations:

.. code-block:: bash

   GenoKit kmer -f small.fa --source-type genome -l 20 \
     -o results/small_kmer -r --engine python

Use Jellyfish for large human inputs. Multiple k values increase preparation
time and disk usage. A catalogue may contain counts internally or refer to
external ``.jf`` payloads; moving only the SQLite file is not always sufficient.

Register an existing dataset
----------------------------

.. code-block:: bash

   GenoKit kmer --register-jf human.fa.kmer.20.jf --source-type genome \
     -o results/kmer -r -j /path/to/jellyfish

``--register-jf`` can be repeated. Registration records existing counts without
recounting the reference. ``--sqlite-db`` selects a catalogue path explicitly.
Register a source and canonical mode consistent with how the dataset was built;
registration is not a conversion of the underlying counts.

Move between platforms
----------------------

Copy the catalogue and its referenced datasets, install Jellyfish for the new
platform, and register each moved dataset using its new path and executable:

.. code-block:: bash

   GenoKit kmer --register-jf /data/kmer/human.fa.kmer.20.jf \
     --source-type genome --sqlite-db /data/kmer/genokit.kmer.sqlite \
     -o /data/kmer -r -j /opt/bio/bin/jellyfish

Use the same source/k combination to update its registration. Re-register all
moved sources/lengths that will be queried. A design command's ``--jellyfish``
can override the stored executable path, but cannot repair a missing dataset.
Retain reference and annotation checksums alongside the catalogue: metadata
about source type alone cannot establish that the assembly/release matches.
