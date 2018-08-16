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
  - spatialindex

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

install spatialindex and rtree

```
-curl -L https://github.com/libspatialindex/libspatialindex/archive/1.8.5.tar.gz | tar xz
cd spatialindex-src-1.8.5
./configure
make
sudo make install
sudo ldconfig
pip install rtree
```
Run the R code:

```
install.packages("argparser")
install.packages("devtools")
library(devtools)
install_bitbucket("chicagoTeam/Chicago", subdir="Chicago")
```
