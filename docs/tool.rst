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


Tools for processing fastq C-HiC files
================================

.. automodule:: tool

   Read truncation
   ===============

   HiCUP truncater
   ---------------
   .. autoclass:: tool.truncater.Truncater
      :members:

   Map and parser reads
   ======================

   fastq2bed
   ---------
   .. autoclass:: tool.fastq2bed.Fastq2bed
      :members:

   Create CHiCAGO input files
   ==========================

   makeRmap
   ---------
   .. autoclass:: tool.makeRmap_Tool.makeRmapFile           
      :members:

   makeBaitmap     
   -----------
   .. autoclass:: tool.makeRmap_Tool.makeRmapFile
      :members:



test_makeRmap_Tool.py')
    params.append('test_makeBaitmap.py')
    params.append('test_makeDesignFiles_Tool.py')
    params.append('test_fastq2bed.py')
    params.append('test_bed2bamchicago.py')
    params.append('test_bam2chicago_Tool.py')
   # params.append('test_runChicago.py')

