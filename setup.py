# -*- coding: utf-8 -*-
import os

def readme():
    with open("README.md") as f:
        long_description = f.read()
        return long_description

from setuptools import setup 
from GenoKit.version import __version__

PACKAGES = [
    "GenoKit",
    "GenoKit.commands",
    "GenoKit.database",
    "GenoKit.utils"
]
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
    packages=PACKAGES,
    license='MIT',
    install_requires = ['argparse>=1.1', 
                        'pandas>=1.2.4', 
                        'gffutils>=0.10.1',
                        'setuptools>=49.2.0',
                        'biopython>=1.78'],
    python_requires=">=3.7.6")
