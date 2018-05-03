# C-HiC

[![Documentation Status](https://readthedocs.org/projects/mg-process-test/badge/?version=latest)](http://mg-process-test.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/mg-process-test.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/mg-process-test) [![Code Health](https://landscape.io/github/Multiscale-Genomics/mg-process-test/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/mg-process-test/master)


This repo contains tools and pipelines for analyzing capture Hi-C data. CHiCAGO algorithm is used for the normalization of the chromatin contacts

# Requirements
- pyenv and pyenv-virtualenv
- Python 2.7.12
- Python Modules:
  - pylint
  - pytest
  - mg-tool-api
  - rpy2
- R >=3.1.2
-R Modules:
  -argparser
  -devtools
  -Chicago
- bedtools
- perl

Installation
------------

Directly from GitHub:

```
cd ${HOME}/code

git clone https://github.com/Multiscale-Genomics/mg-process-test.git

cd mg-process-test
```

Create the Python environment

```
pyenv-virtualenv 2.7.12 mg-process-test
pyenv activate mg-process-test
pip install -e .
pip install -r requirements.txt
```
Run the R code:

```
install.packages("argparser")
install.packages("devtools")
library(devtools)
install_bitbucket("chicagoTeam/Chicago", subdir="Chicago")
```
