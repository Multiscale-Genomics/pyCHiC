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

Requirements and Installation
=============================

Requirements
------------

Software
^^^^^^^^

- Python 2.7.12+
- R >=3.1.2
- bedtools
- perl


Python Modules
^^^^^^^^^^^^^^

- mg-tool-api
- pylint
- pytest

R Modules
^^^^^^^^^
- argparser
- devtools
- Chicago

To Run runChicago.py and process_runChicago.py, the R script runChicago.R from  https://bitbucket.org/chicagoTeam/chicago/src/ceffddda8ea392a1e84e4db9593f8fc35ac88048/chicagoTools/?at=master
should be downloded and added to PATH.

Installation
------------
Directly from GitHub:

.. code-block:: none
   :linenos:

   git clone https://github.com/pabloacera/C-HiC.git

Using pip:

.. code-block:: none
   :linenos:

   pip install git+https://github.com/pabloacera/C-HiC.git

Install R modules, use the following R code:

  install.packages("argparser")
  install.packages("devtools") 
  library(devtools)
  install_bitbucket("chicagoTeam/Chicago", subdir="Chicago")
	
