# -*- coding: utf-8 -*-
from pathlib import Path

from setuptools import find_packages, setup

from GenoKit.version import __version__


ROOT = Path(__file__).resolve().parent


def readme():
    return (ROOT / "README.md").read_text(encoding="utf-8")


setup(
    name='GenoKit',
    version=__version__,
    keywords='genomic feature, extract',
    description='A versatile command-line toolkit for comprehensive genomic feature extraction, analysis, and visualization',
    long_description=readme(),
    long_description_content_type='text/markdown',
    entry_points = {'console_scripts': [
                       'GenoKit=GenoKit.command_genokit:main',
                       'GenoKitGB=GenoKit.command_genokitgb:main'
                   ]},
    author='zhusitao',
    author_email='zhusitao1990@163.com',
    url='https://github.com/SitaoZ/GenoKit.git',
    include_package_data=True, # done via MANIFEST.in under setuptools
    package_data={'GenoKit.example': ['*.csv']},
    packages=find_packages(exclude=("test", "test.*")),
    license='MIT',
    install_requires = ['pandas>=2.2.3',
                        'setuptools>=72.1.0',
                        'biopython>=1.86',
                        'python-louvain>=0.16',
                        'python-circos>=0.3.0',
                        'tabulate>=0.9.0',
                        'tqdm>=4.0',
                        'Django>=5.2.8',
                        'pyfaidx>=0.9',
                        'matplotlib>=3.10.0',
                        'requests>=2.32.3',
                        'networkx>=3.6.1'],
    python_requires=">=3.10")
