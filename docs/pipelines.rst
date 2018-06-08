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

Pipelines
=========

Truncate reads
--------------
.. automodule::


   Running from the command line
   =============================

   Parameters
   ----------
   config : file
      Location of the config file for the workflow
   in_metadata : file
      Location of the input list of files required by the process
   out_metadata : file
      Location of the output results.json file for returned files

   Returns
   -------
   output : file
      Text file with a single entry

   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):

   .. code-block:: none
      :linenos:

      cd /home/compss/code/mg-process-test
      runcompss --lang=python process_test.py --config /home/compss/code/mg-process-test/tool_config/process_test.json --in_metadata /home/compss/code/mg-process-test/tests/json/input_test.json --out_metadata /home/compss/code/mg-process-test/tests/results.json

   Methods
   =======
   .. autoclass:: process_test.process_test
      :members:






















Map and parse C-HiC reads
-------------------------
.. automodule:: process_fastq2bed

   This pipeline will take as input two fastq files, RE sites, the genome indexed with GEM and the same genome
   in FASTA file. This pipeline uses TADbit to map, filter and produce a bed file that will be used later on to
   produce bam file compatible with CHiCAGO algorithm. More information about filtering and mapping https://3dgenomes.github.io/TADbit/

   Running from the command line
   =============================

   Parameters
   ----------
   config : str
      Configuration JSON file
   in_metadata : str
      Location of input JSON metadata for files
   out_metadata : str
      Location of output JSON metadata for files

   Returns
   -------
   Wd : folders and files
      path to the working directory
      where the output files are


   Example
   -------
   REQUIREMENT - Needs two fastq files single end, FASTA genome and GEM indexed genome

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_fastq2bed.py \
         --config tests/json/config_fastq2bed.json \
         --in_metadata tests/json/input_biobambam.json \
         --out_metadata tests/json/output_biobambam.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_biobambam.py         \
            --config tests/json/config_biobambam.json \
            --in_metadata tests/json/input_biobambam.json \
            --out_metadata tests/json/output_biobambam.json


   Methods
   =======
   .. autoclass:: process_biobambam.process_biobambam
      :members:






















process_baitmap.py
-rw-rw-r--  1 pablo pablo  4861 May 21 12:51 process_baitmap.pyc
-rw-rw-r--  1 pablo pablo  4735 May 22 10:54 process_bam2chicago_Tool.py
-rw-rw-r--  1 pablo pablo  4591 May 21 17:23 process_bam2chicago_Tool.pyc
-rw-rw-r--  1 pablo pablo  4489 May 22 11:07 process_bed2bam.py
-rw-rw-r--  1 pablo pablo  4562 May 21 16:41 process_bed2bam.pyc
-rw-rw-r--  1 pablo pablo  6270 May 24 15:05 process_chicago_CHiC.py
-rw-rw-r--  1 pablo pablo  5714 May 23 09:45 process_chicago_CHiC.pyc
-rw-rw-r--  1 pablo pablo  4805 May 22 15:29 process_Design.py
-rw-rw-r--  1 pablo pablo  4791 May 21 12:34 process_Design.pyc
-rw-rw-r--  1 pablo pablo  5345 May 22 15:20 process_fastq2bed.py
-rw-rw-r--  1 pablo pablo  5251 May 21 15:52 process_fastq2bed.pyc
-rw-rw-r--  1 pablo pablo  4915 May 22 15:52 process_rmap.py
-rw-rw-r--  1 pablo pablo  4810 May 21 12:36 process_rmap.pyc
-rw-rw-r--  1 pablo pablo  4832 May 22 15:20 process_runChicago.py
-rw-rw-r--  1 pablo pablo  4716 May 21 12:10 process_runChicago.pyc

