.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Full Installation
=================

The following document is for the full installation of all software required by
the C-HiC module and all programmes that it uses. The document has
been written with Ubuntu Linux, although many of the commands are cross
platform (\*nix) complient.

If you already have certain packages installed feel free to skip over certain
steps. Likewise the bin, lib and code directories are relative to the home dir;
if this is not the case for your system then make the required changes when
running these commands.

Setup the System Environment
----------------------------

.. code-block:: none
   :linenos:

   sudo apt-get install -y make build-essential libssl-dev zlib1g-dev       \\
   libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \\
   libncursesw5-dev xz-utils tk-dev unzip mcl libgtk2.0-dev r-base-core     \\
   libcurl4-gnutls-dev python-rpy2 git libtbb2 pigz liblzma-dev libhdf5-dev \\
   texlive-latex-base

   cd ${HOME}
   mkdir bin lib code
   echo 'export PATH="${HOME}/bin:$PATH"' >> ~/.bash_profile

Setup pyenv and pyenv-virtualenv
--------------------------------

This is required for managing the version of Python and the installation
environment for the Python modules so that they can be installed in the user
space.

.. code-block:: none
   :linenos:

   git clone https://github.com/pyenv/pyenv.git ~/.pyenv
   echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
   echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
   echo 'eval "$(pyenv init -)"' >> ~/.bash_profile

   # Add the .bash_profile to your .bashrc file
   echo 'source ~/.bash_profile"' >> ~/.bashrc

   git clone https://github.com/pyenv/pyenv-virtualenv.git ${PYENV_ROOT}/plugins/pyenv-virtualenv

   pyenv install 2.7.12
   pyenv virtualenv 2.7.12 C-HiC

   # Python 3 environment required by iNPS
   pyenv install 3.5.3
   ln -s ${HOME}/.pyenv/versions/3.5.3/bin/python ${HOME}/bin/py3

Installation Process
--------------------


bedtools
^^^^^^^^
.. code-block:: none
   :linenos:

   sudo apt-get install bedtools

Bowtie2 Aligner
^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   wget --max-redirect 1 https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip
   unzip bowtie2-2.3.4-linux-x86_64.zip

HiCUP
^^^^^
.. code-block:: none
   :linenos:
    
    cd ${HOME}/lib
    wget https://www.bioinformatics.babraham.ac.uk/projects/hicup/hicup_v0.6.1.tar.gz
    tar -xzf hicup_v0.6.1.tar.gz
    cd hicup_v0.6.1
    chmod a+x *
    #modify hicup_truncater for delete the time of the summary file



BWA Sequence Aligner
^^^^^^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/lh3/bwa.git
   cd bwa
   make

SAMtools
^^^^^^^^

.. code-block:: none
   :linenos:

   cd ${HOME}/lib
   git clone https://github.com/samtools/htslib.git
   cd htslib
   autoheader
   autoconf
   ./configure --prefix=${HOME}/lib/htslib
   make
   make install

   cd ${HOME}/lib
   git clone https://github.com/samtools/samtools.git
   cd samtools
   autoheader
   autoconf -Wno-syntax
   ./configure --prefix=${HOME}/lib/samtools
   make
   make install

Install CHiCAGO
^^^^^^^^^^^^^^^

.. code-block:: none
   :linenos:

   sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
   sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
   sudo apt-get update -qq
   sudo apt-get install r-base-core
   sudo apt-get install python-rpy2


   cd ${HOME}/lib
   sudo apt-get install libtbb-dev
   sudo apt-get install libssl-dev
   cd ${HOME}/C-HiC/
   echo "R_LIB=${HOME}/R" > ${HOME}/.Renviron
   echo "options(repos = c(CRAN = 'http://mirrors.ebi.ac.uk/CRAN/'))" > ${HOME}/.Rprofile
   echo ".libPaths('~/R')" >> ${HOME}/.Rprofile
   echo 'message("Using library:", .libPaths()[1])' >> ${HOME}/.Rprofile
   sudo Rscript scripts/install_packages.R


Setup the symlinks
------------------

.. code-block:: none
   :linenos:

   cd ${HOME}/bin
    
      

   ln -s ${HOME}/lib/hicup_v0.6.1/* ${HOME}/bin/

   ln -s ${HOME}/lib/bwa/bwa bwa

   ln -s ${HOME}/lib/htslib/bin/bgzip bgzip
   ln -s ${HOME}/lib/htslib/bin/htsfile htsfile
   ln -s ${HOME}/lib/htslib/bin/tabix tabix

   ln -s ${HOME}/lib/samtools/bin/ace2sam ace2sam
   ln -s ${HOME}/lib/samtools/bin/blast2sam.pl blast2sam.pl
   ln -s ${HOME}/lib/samtools/bin/bowtie2sam.pl bowtie2sam.pl
   ln -s ${HOME}/lib/samtools/bin/export2sam.pl export2sam.pl
   ln -s ${HOME}/lib/samtools/bin/interpolate_sam.pl interpolate_sam.pl
   ln -s ${HOME}/lib/samtools/bin/maq2sam-long maq2sam-long
   ln -s ${HOME}/lib/samtools/bin/maq2sam-short maq2sam-short
   ln -s ${HOME}/lib/samtools/bin/md5fa md5fa
   ln -s ${HOME}/lib/samtools/bin/md5sum-lite md5sum-lite
   ln -s ${HOME}/lib/samtools/bin/novo2sam.pl novo2sam.pl
   ln -s ${HOME}/lib/samtools/bin/plot-bamstats plot-bamstats
   ln -s ${HOME}/lib/samtools/bin/psl2sam.pl psl2sam.pl
   ln -s ${HOME}/lib/samtools/bin/sam2vcf.pl sam2vcf.pl
   ln -s ${HOME}/lib/samtools/bin/samtools samtools
   ln -s ${HOME}/lib/samtools/bin/samtools.pl samtools.pl
   ln -s ${HOME}/lib/samtools/bin/seq_cache_populate.pl seq_cache_populate.pl
   ln -s ${HOME}/lib/samtools/bin/soap2sam.pl soap2sam.pl
   ln -s ${HOME}/lib/samtools/bin/varfilter.py varfilter.py
   ln -s ${HOME}/lib/samtools/bin/wgsim wgsim
   ln -s ${HOME}/lib/samtools/bin/wgsim_eval.pl wgsim_eval.pl
   ln -s ${HOME}/lib/samtools/bin/zoom2sam.pl zoom2sam.pl

Prepare the Python Environment
------------------------------

Install APIs and Pipelines
^^^^^^^^^^^^^^^^^^^^^^^^^^

Checkout the code for the DM API and the C-HiC pipelines:

.. code-block:: none
   :linenos:

   cd ${HOME}/code
   pyenv activate C-HiC
   pip install git+https://github.com/Multiscale-Genomics/mg-dm-api.git
   pip install git+https://github.com/Multiscale-Genomics/mg-tool-api.git
   pip install git+https://github.com/Multiscale-Genomics/mg-process-fastq.git


   git clone https://github.com/pabloacera/C-HiC.git
   cd C-HiC
   pip install -e .
   pip install -r requirements.txt
   pip install dill
