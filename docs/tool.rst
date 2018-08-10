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

   Map and parser reads
   ======================

   hicup_tool
   ---------
   .. autoclass:: CHiC.tool.hicup_tool.hicup
      :members:

   Create CHiCAGO input files
   ==========================

   makeRmap
   ---------
   .. autoclass:: CHiC.tool.makeRmap_Tool.makeRmapFile
      :members:

   makeBaitmap
   -----------
   .. autoclass:: CHiC.tool.makeRmap_Tool.makeRmapFile
      :members:

   makeDesignFiles
   ---------------
   .. autoclass:: CHiC.tool.makeDesignFiles_Tool.makeDesignFilesTool
      :members:

   Convert bam files into chicago input
   =====================================

   bam2chicago
   -----------
   .. autoclass:: CHiC.tool.bam2chicago_Tool.bam2chicagoTool
      :members:

   Normalize data and call C-HiC peaks
   ===================================

   CHiCAGO
   -------
   .. autoclass:: CHiC.tool.runChicago.ChicagoTool
      :members: