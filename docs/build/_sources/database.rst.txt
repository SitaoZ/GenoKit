Annotation databases and utilities
==================================

Create and reuse a database
---------------------------

.. code-block:: bash

   GenoKit create -s gtf -a human.gtf -g human.fa -o results -p human
   GenoKit create -s gff -a annotation.gff3.gz -g genome.fa -o results -p sample

The default products are ``results/human.gtf.sqlite`` and
``results/sample.gff.sqlite``. The database stores annotation structure and
attributes; keep the reference FASTA separately. The default builder inserts
records incrementally and the reader loads gene models lazily. The implementation
uses GenoKit's SQLite schema, not a gffutils database.

The builder checks annotation syntax against the chosen style. Missing exon
structure, inconsistent parent IDs or a reference with different contig names
should be addressed before downstream design.

Legacy backend
--------------

.. code-block:: bash

   GenoKit create -s gtf -a human.gtf -g human.fa -o results -p human \
     --database-format pickle

This writes ``human.gtf.pkl``. Readers detect SQLite by its file signature and
otherwise support historical pickle representations, including older ``.bin``
files. Pickle deserializes the full stored object and is less suitable for large
annotations; use only trusted pickle inputs. Renaming another SQLite schema
will not convert it into a GenoKit annotation database.

Index handling
--------------

``create`` normally checks/builds the FASTA index. ``--skip-fasta-index`` skips
that preparation, not the indexing requirement of later sequence access.
When moving to another machine, move the annotation database and matching
reference, and give their new paths with ``-d`` and ``-g``. Annotation databases
and k-mer catalogues have different schemas and are not interchangeable.

Statistics and record selection
-------------------------------

.. code-block:: bash

   GenoKit stat -d results/human.gtf.sqlite -g human.fa -s gtf -o results/stat.csv
   GenoKit fasta -f human.fa -i chromosomes.txt -o results/chromosomes.fa
   GenoKit iupac -n
   GenoKit iupac -a
   GenoKit iupac -c

``chromosomes.txt`` contains one FASTA record ID per line, with no header.
``fasta`` selects complete records rather than an interval such as chr1:100-200.
The three ``iupac`` examples print nucleotide codes, amino-acid codes and the
codon table respectively. ``stat`` can print statistics when ``-o`` is omitted.
