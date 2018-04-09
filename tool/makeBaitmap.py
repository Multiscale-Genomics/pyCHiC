"""
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
"""

from __future__ import print_funtion

import os
import subprocess
import sys

from utils import logger


try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

##################################################

class makeBaitmapTool(Tool):
   """
   This tool use probe capture sequences as input.
   Then using bwa can tell which baits correspond
   to the probes
   """
   def __init__(self, configuration = None):
      """
      Initialise the tool with its configuration
      Parameters:
      -----------
      configuration: dict
         parameters to run the tool
      """

      logger.info("initialising makeBaitmapTool")
      Tool.__init__(self)

      if configuration is None:
         configuration = {}

      self.configuration.update(configuration)


   def bwa_for_probes(self, genome_index, probes_fa, out):
      """
      This function run bwa using an index genome and a probes file
      in fasta format. bwa is used as single end and with high
      gap penalty and missmacht score
      Parameters:
      -----------
      genome_index: str
         path to the indexed genome. This genome should
         be the same used to generate the .rmap file
      probes_fa: str
         path to probes files in fasta format,
         every sequences representing one
      out: str
         Name of the output file
      Return:
      ------
         bool
      """

      args = ["bwa mem -O 100 -B 20", genome_index,
      probes_fa, out]

      logger.info("bwa_for_probes CMD: " " ".join(args))

      process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr= subprocess.PIPE)
      process.wait()
      proc_out, proc_err = process.communicate()

      if os.path.getsize(out + ".sam") is True:
         return True
      else:
         logger.fatal("bwa failed to generate sam file")
         logger.fatal("bwa stdout" + proc_out)
         logger.fatal("bwa stderr" + proc_err)
         return False

   def reformat_to_baitmap(self, sam_file, rmap, baitmap_file):
      """
      This function take the sam file, output of bwa
      and the rmap file, and output a baitmap file
      Parameters:
      -----------
      sam_file : str
         path to output file from bwa_for_probes
      rmap: str
         complete path to .rmap file
      baitmap_file: str
         complete path and name of the baitmap file
      """

      with open("rmap", "r") as file:
         for line in file:
            line = line.strip().split("\t")




      with open("sam_file", "r") as file_in:
         with open("baitmap_file", "a") as file_put:
            for line in file_in:
               line = line.rstrip().split("\t")
               if line[0][0] != "@":
