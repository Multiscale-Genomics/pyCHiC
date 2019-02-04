# pyCHiC

[![Documentation Status](https://readthedocs.org/projects/capture-chi-c/badge/?version=latest)](https://capture-chi-c.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/pyCHiC.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/pyCHiC)


This repository contains pipelines for analyzing capture Hi-C data. CHiCAGO algorithm is used for the normalization of chromatin contacts

# Requirements
- pyenv and pyenv-virtualenv
- Python 3.6
- Python Modules:
  - pylint
  - pytest
  - mg-tool-api
  - rpy2
  - matplotlib==2.2.3
  - pandas==0.23.4
  - rtree
  - rpy2
  - numpy==1.14.5
  - pandas
  - scipy==1.1.0
  - multiprocess
  - matplotlib

- R >=3.1.2
-R Modules:
  -argparser
  -survival
  -Hmisc
  -devtools
  -Chicago
- bedtools
- perl
- spatialindex
- bowtie2
- hicup


Installation
------------

For a guide to the full installation procedure the see [ReadTheDocs](https://capture-chi-c.readthedocs.io/en/latest/).

Directly from GitHub:

.. code-block:: none
   :linenos:

   cd ${HOME}/code

   git clone https://github.com/Multiscale-Genomics/pyCHiC.git

   cd CHi-C

Create the Python environment

.. code-block:: none
   :linenos:

   pyenv-virtualenv 2.7.10 CHi-C
   pip install --editable .

Tests
-----

Tests must be runned from the top directory.
Test for single tools are runned with pytest, example: $pytest tests/test_rmap_tool.py
There is an order to run single tests:
  - test_rmap_tool.py
  - test_baitmap.py
  - test_design.py
  - test_hicup.py
  - test_bam2chicago.py
  - test_pyCHiC.py

To run all test type: $python tests/test_toolchains.py

IMPORTANT! From the top directory type ./tidy_data before running tests/test_toolchains.py

To test run_pyCHiC_pipeline:
python run_pyCHiC_pipeline.py --bowtie_idx tests/data/test_baitmap/chr21_hg19.fa.bt2.tar.gz  --genome_fa tests/data/test_baitmap/chr21_hg19.fa --probes_fa tests/data/test_baitmap/h19_promoter.fa --fastq1 tests/data/test_truncater/SRR3535023_1_chr21_new.fastq --fastq2 tests/data/test_truncater/SRR3535023_2_chr21_new.fastq --execution . --RE_name HindIII --RE_sequence A^AGCTT --execution test_run --minNPerBait 1 --minProxOEPerBin 1 --minProxB2BPerBin 1 --techNoise_minBaitsPerBin 1 --cutoff 1
