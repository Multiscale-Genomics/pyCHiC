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

Architectural Design Record
===========================

25-09-2018 handling_chr_header branch merge with master
-------------------------------------------------------

This rmap_tool.py from this branch take the chromosome format from the used the reference genome and
output a file with two columns, dictionary like with number of the chromosome and the name of the chromsome from the reference genome. example
1 chr1
2 chr2
3 chr3
ect...

This file is passed to the makeBaitmap.py script and generate the .batimap file with the corresponding chromsome name. This is necesary as the rtrees used in makeBaitmap.py needs an integer instead of "chr" or any other format.

15-10-2018 mm_mods_for_makebaitmaps branch merge with master
------------------------------------------------------------

This branch contains some modifications from Mark to solve issues with pyCOMPSs regarding makeBaitmap.py tool


25-1-2019 Creating branch pyCOMPSs
----------------------------------
This is the branch containing pyCOMPSs decorators, and is working in COMPSs. This Branch was not
merged with master because the code contain modifications that will made the code run slower in local
(mostly the function getEtaBar is not running in parallel due to issues in COMPSs)


25-2-2019 Creating branch called folder_name
--------------------------------------------

This branch rename the folder CHiC to pyCHiC. The path from the pyCHiC master and pyCOMPSs branches to acces the tools was pyCHiC/CHiC/tool/.. This branch rename CHiC and convert the path to pyCHiC/pyCHiC/tool/...


