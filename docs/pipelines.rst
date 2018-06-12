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
.. automodule:: process_truncater

   This tool check single end Capture-HiC fastq reads and eliminates bases pairs that are
   beyond the RE site.

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
   output : files
     Truncated fastq files.

   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_truncater.py         \
            --config tests/json/config_truncater.json \
            --in_metadata tests/json/input_truncater.json \
            --out_metadata tests/json/output_truncater.json

   Methods
   =======
   .. autoclass:: process_truncater.process_truncater
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
   REQUIREMENT - Needs two fastq files single end, FASTA genome and GEM indexed genome.

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_fastq2bed.py \
         --config tests/json/config_fastq2bed.json \
         --in_metadata tests/json/input_fastq2bed.json \
         --out_metadata tests/json/output_fastq2bed.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_fastq2bed.py         \
            --config tests/json/config_fastq2bed.json \
            --in_metadata tests/json/input_fastq2bed.json \
            --out_metadata tests/json/output_fastq2bed.json

   Methods
   =======
   .. autoclass:: process_fastq2bed.process_Fastq2bed
      :members:


Create CHiCAGO input RMAP
-------------------------
.. automodule:: process_rmap

   This pipeline creates the .rmap file, one of the inputs of CHiCAGO. The file Consisting on
    <chr> <start> <end> <numeric ID>. Is a virtual digest of the genome using a RE.

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
   output_files : .rmap file
   Rtree_files: rtree file with information
      about the RE fragments in the genome. It is
      used for the process_baitmap.py

   Example
   -------
   REQUIREMENT - Needs FASTA file of the gneome, and a RE in the config file.

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_rmap.py \
         --config tests/json/config_rmap.json \
         --in_metadata tests/json/input_rmap.json \
         --out_metadata tests/json/output_rmap.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_rmap.py         \
            --config tests/json/config_fastq2bed.json \
            --in_metadata tests/json/input_fastq2bed.json \
            --out_metadata tests/json/output_fastq2bed.json

   Methods
   =======
   .. autoclass:: process_fastq2bed.process_Fastq2bed
      :members:



Create CHiCAGO input BAITMAP
-------------------------
.. automodule:: process_baitmap

   This pipeline creates the .baitmap file, one of the inputs of CHiCAGO. The file Consisting on
   <chr> <start> <end> <numeric ID> <annotation> Is a subset of the rmap file, containing the RE fragments that
   overlap with baits provided by the user. Baits are the RE fragments that are capture during the experimental protocol.

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
   out_baitmap : .baitmap file
   out_sam : .sam file
      used to generate .baitmap

   Example
   -------
   REQUIREMENT - Needs rtree file generated using process_rmap.py,
                 genome indexed using bwa,
                 File with the used probes,

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_baitmap.py \
         --config tests/json/config_baitmap.json \
         --in_metadata tests/json/input_baitmap.json \
         --out_metadata tests/json/output_baitmap.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_rmap.py         \
            --config tests/json/config_baitmap.json \
            --in_metadata tests/json/input_baitmap.json \
            --out_metadata tests/json/output_baitmap.json

   Methods
   =======
   .. autoclass:: process_baitmap.generate_CHiCAGO_baitmap
      :members:



Create CHiCAGO input Design files
---------------------------------
.. automodule:: process_Desgin

   This script use as input .rmap and .baitmap files and generate the Design files.

   NPerBin file (.npb): <baitID> <Total no. valid restriction fragments in distance bin 1> ... <Total no. valid
   restriction fragments in distance bin N>,
   where the bins map within the "proximal" distance range from each bait (0; maxLBrownEst] and bin size is defined by
   the binsize parameter.
   NBaitsPerBin file (.nbpb): <otherEndID> <Total no. valid baits in distance bin 1> ... <Total no. valid baits in
   distance bin N>,
   where the bins map within the "proximal" distance range from each other end (0; maxLBrownEst] and bin size is defined
   by the binsize parameter.
   Proximal Other End (ProxOE) file (.poe): <baitID> <otherEndID> <absolute distance>
   for all combinations of baits and other ends that map within the "proximal" distance range from each other (0;
   maxLBrownEst].
   Data in each file is preceded by a comment line listing the input parameters used to generate them.


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
   designDir : dir with the design output files
   out_sam : .sam file
      used to generate .baitmap

   Example
   -------
   REQUIREMENT - Needs .rmap and .baitmap files.

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_Design.py \
         --config tests/json/config_Design.json \
         --in_metadata tests/json/input_Design.json \
         --out_metadata tests/json/output_Design.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_Design.py         \
            --config tests/json/config_Design.json \
            --in_metadata tests/json/input_Design.json \
            --out_metadata tests/json/output_Design.json

   Methods
   =======
   .. autoclass:: process_Design.process_makeDesign
      :members:



Convert bed files into BAM chicago-friendly files
-------------------------------------------------
.. automodule:: process_bed2bam

   This pipeline convert the output of process_fastq2bed.py file called
   valid_r1-r2_intersection_b51cdf1282.tsv to a bam file compatible with
   bam2chicago.py

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
   bam_out: bam output file

   Example
   -------
   REQUIREMENT - Needs bed file produced by process_fastq2bed

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_bed2bam.py \
         --config tests/json/config_bed2bam.json \
         --in_metadata tests/json/input_bed2bam.json \
         --out_metadata tests/json/output_bed2bam.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/
   software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_Design.py         \
            --config tests/json/config_bed2bam.json \
            --in_metadata tests/json/input_bed2bam.json \
            --out_metadata tests/json/output_bed2bam.json

   Methods
   =======
   .. autoclass:: process_bed2bam.process_bed2bam
      :members:



Convert BAM file into chicago input files .chinput
---------------------------------------------------
.. automodule:: process_bam2chicago_Tool

   This pipeline convert the output of process_bed2bam.py BAM file to
   a .chinput file, input for process_runChicago.py

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
   chrRMAP : .rmap file with chr# format
   chrBAITMAP : .baitmap file with chr# format
   sample_name : .chinput output

   Example
   -------
   REQUIREMENT - Needs BAM file produced by process_bed2bam.py
                 Needs a .rmap file
                 Needs a .baitmap file

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_bam2chicago_Tool.py \
         --config tests/json/config_bam2chicago_Tool.json \
         --in_metadata tests/json/input_bam2chicago_Tool.json \
         --out_metadata tests/json/output_bam2chicago_Tool.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/
   software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_bam2chicago_Tool.py         \
            --config tests/json/config_bam2chicago_Tool.json \
            --in_metadata tests/json/input_bam2chicago_Tool.json \
            --out_metadata tests/json/output_bam2chicago_Tool.json

   Methods
   =======
   .. autoclass:: process_bam2chicago_Tool.process_bam2chicago
      :members:



Data normalization and peak calling
-----------------------------------
.. automodule:: process_runChicago

   This pipeline runs the normalization of the data and call the real
   chomatine interactions

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
   output_dir: directory with all output folders and files

   Example
   -------
   REQUIREMENT - Needs at least one .chinput file
                 Needs a config file with:
                     - settings file
                     - design dir:
                           .rmap
                           .baitmap
                           .npb
                           .nbpb
                           .poe

   When running the pipeline on a local machine without COMPSs:

   .. code-block:: none
      :linenos:

      python process_runChicago.py \
         --config tests/json/config_runChicago.json \
         --in_metadata tests/json/input_runChicago.json \
         --out_metadata tests/json/output_runChicago.json \
         --local

   When using a local version of the [COMPS virtual machine](https://www.bsc.es/research-and-development/
   software-and-apps/software-list/comp-superscalar/):

   .. code-block:: none
      :linenos:

      runcompss                     \
         --lang=python              \
         --library_path=${HOME}/bin \
         --pythonpath=/<pyenv_virtenv_dir>/lib/python2.7/site-packages/ \
         --log_level=debug          \
         process_runChicago.py         \
            --config tests/json/config_runChicago.json \
            --in_metadata tests/json/input_runChicago.json \
            --out_metadata tests/json/output_runChicago.json

   Methods
   =======
   .. autoclass:: process_runChicago.process_runChicago
      :members: