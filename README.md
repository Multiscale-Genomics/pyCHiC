# CHi-C

[![Documentation Status](https://readthedocs.org/projects/CHi-C/badge/?version=latest)](http://CHi-C.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/CHi-C.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/CHi-C) [![Code Health](https://landscape.io/github/Multiscale-Genomics/CHi-C/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/CHi-C/master)


This repository contains pipelines for analyzing capture Hi-C data. CHiCAGO algorithm is used for the normalization of chromatin contacts

# Requirements
- pyenv and pyenv-virtualenv
- Python 2.7.12
- Python Modules:
  - pylint
  - pytest
  - mg-tool-api
  - rpy2
  - matplotlib
  - pandas
  - rtree


- R >=3.1.2
-R Modules:
  -argparser
  -devtools
  -Chicago
- bedtools
- perl
- spatialindex
- bowtie2
- hicup
- BWA


Installation
------------

For a guide to the full installation procedure the see [ReadTheDocs](http://CHi-C.readthedocs.io).

Directly from GitHub:

.. code-block:: none
   :linenos:

   cd ${HOME}/code

   git clone https://github.com/Multiscale-Genomics/CHi-C.git

   cd CHi-C

Create the Python environment

.. code-block:: none
   :linenos:

   pyenv-virtualenv 2.7.10 CHi-C
   pip install --editable .
