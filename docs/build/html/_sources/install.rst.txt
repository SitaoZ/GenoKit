Installation
============

Choose an isolated Python environment
-------------------------------------

The source metadata declares Python >=3.10, but dependency requirements also
apply. In particular, the declared NetworkX >=3.6.1 requirement makes Python
3.11 or newer the practical choice for this dependency set. Python 3.8 guidance
from older GenoKit documentation does not describe this source release.

.. code-block:: bash

   python3.11 -m venv genokit-env
   source genokit-env/bin/activate
   python -m pip install --upgrade pip
   python -m pip install GenoKit
   GenoKit --version
   GenoKit -h

On Windows, activate a venv with ``genokit-env\Scripts\activate`` in cmd or
``genokit-env\Scripts\Activate.ps1`` in PowerShell. Examples elsewhere use
Bash syntax; adapt line continuations and executable paths to your shell.

Install the current source
--------------------------

.. code-block:: bash

   git clone https://github.com/SitaoZ/GenoKit.git
   cd GenoKit
   python -m pip install .

Use ``python -m pip install -e .`` for an editable checkout. Prefer pip over
``python setup.py install``. The PyPI release and a development checkout can
have different features; record which one was used.

Declared Python dependencies
----------------------------

These are the minimum versions in the current ``setup.py``; they are not a
complete lockfile for a benchmark environment.

.. list-table::
   :header-rows: 1
   :widths: 45 25

   * - Package
     - Minimum version
   * - pandas
     - 2.2.3
   * - setuptools
     - 72.1.0
   * - biopython
     - 1.86
   * - python-louvain
     - 0.16
   * - python-circos
     - 0.3.0
   * - tabulate
     - 0.9.0
   * - tqdm
     - 4.0
   * - Django
     - 5.2.8
   * - pyfaidx
     - 0.9
   * - matplotlib
     - 3.10.0
   * - requests
     - 2.32.3
   * - networkx
     - 3.6.1

gffutils is not a declared dependency of the current annotation backend.
R, GenomicFeatures, AGAT and BEDTools are comparison tools, not prerequisites
for ordinary GenoKit extraction.

Optional external tools
-----------------------

Jellyfish is required when building/querying Jellyfish-backed k-mer datasets.
Install it for the destination operating system, then provide its executable
with ``GenoKit kmer -j /path/to/jellyfish`` or a design command's
``--jellyfish /path/to/jellyfish``. The Python k-mer engine is an alternative
for small inputs; see :doc:`kmer` for its limits.

For compressed reference FASTA, use BGZF rather than ordinary gzip. A bgzip
executable can prepare that input, but plain FASTA also works. NCBI access is
needed for HGVS reference downloads unless a matching local record is supplied.
Sphinx is needed only to build the documentation, not to run GenoKit.

Check the installation
----------------------

.. code-block:: bash

   python -m pip check
   python -m pip show GenoKit
   GenoKit create -h
   GenoKit primer -h

If the command resolves to another environment, inspect ``command -v GenoKit``
and ``python -m pip --version``. For an editable installation, the checkout must
remain available. See :doc:`faq` for environment and plotting issues.
