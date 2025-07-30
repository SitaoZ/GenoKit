Installation
============

git
------------------
Using `git`::

    git clone https://github.com/SitaoZ/GenoKit.git
    cd GenoKit
    python setup.py install

pip 
---

Or `pip`::

    pip install GenoKit


`GenoKit` is `tested <https://github.com/SitaoZ/GenoKit>`_ with Python 3.8

Optional requirements
---------------------

* python >= 3.7.6 `python <https://www.python.org/>`_
* pandas >= 1.2.4 `pandas <https://pandas.pydata.org/docs/>`_
* gffutils >= 0.10.1 `gffutils <https://pythonhosted.org/gffutils/>`_
* setuptools >= 49.2.0 `setuptools <https://pypi.org/project/setuptools/>`_
* BioPython >= 1.78 `biopython <https://biopython.org/wiki/Documentation/>`_

Install them all with `conda`::

    conda install --channel conda-forge --channel python pandas gffutils setuptools biopython
