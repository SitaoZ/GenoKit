Reproducible workflows and performance
======================================

Record the inputs and environment
---------------------------------

Keep the reference assembly, annotation provider/release, file checksums,
GenoKit source revision or release, commands, dependency versions and output
schemas together. A development source version can differ from installed
package metadata, so include both command output and package information.

.. code-block:: bash

   GenoKit --version > genokit_version.txt
   python --version > python_version.txt
   python -m pip freeze > requirements_run.txt
   python -m pip show GenoKit > genokit_package.txt

For Linux, ``sha256sum genome.fa annotation.gtf`` records input checksums;
on macOS, use ``shasum -a 256``. Archive a copy of small configuration scripts
rather than only a command history containing machine-specific paths.

Design references
-----------------

For k-mer screening, record source type, exact k values, canonical counting,
engine, reference/annotation release and Jellyfish version. For cross-platform
work, rebuild executable environments for the destination OS; do not copy Mac
Python/Perl binaries into Linux. Move data and register dataset/executable paths
as described in :doc:`kmer`.

Measure equivalent tasks
------------------------

Choose whether a measurement covers preparation plus extraction or extraction
from an existing database. Building one program's index inside a timed run but
reusing another's is not an equivalent comparison. Define splicing, strand,
stop-codon inclusion, duplicate records, flanking-window sizes and chromosome
boundary behavior before comparing sequences or performance.

.. code-block:: bash

   python -m pip install memory-profiler psutil
   mprof run --include-children --output transcript_memory.dat \
     GenoKit transcript -d results/human.gtf.sqlite -g human.fa -s gtf \
     -f fasta -o results/transcript_profiled.fa

Run replicates in distinct directories. Record hardware, OS, thread settings,
filesystem/cache conditions, output format and profiler sampling interval.
Sampled memory is not necessarily the instantaneous peak; explain the metric
used in a manuscript. Keep profiling tools in the environment record too.

Do not infer biological correctness from runtime alone. Compare sequence IDs,
lengths and normalized sequence content using equivalent feature definitions.
Empty output, duplicated transcript-associated exons and differing CDS stop
rules can otherwise make a fast run appear incorrectly favorable.
