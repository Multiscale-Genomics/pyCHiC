
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

from __future__ import print_function

import sys
import os
import subprocess
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

class bed2bam(Tool):
    """
    This class contain functions to convert a bed file output from 
    fatq2bed.py to bam file compatible with CHiCAGO
    """
    def __init__(self, configuration=None):
        """
        Initiate the tool

        Parameters
        ----------
        configuration: dict
         contain info to run the functions of the class
        """
        logger.info("Initiating bed2chicago")

    def wrapper_bed2bam(self, bed, ncpus, output):
        """
        This function runs the script from_bed_to_BAM_for_chicago.py
        tha convert a bed file to a bam file compatible with CHiCAGO
        Samtools is necesary for the sccript to run

        Parameters
        ----------
        bed: str
            path to the bed file
        ncpus: str
            Number of cpus to run the script
        output: str
            path to the output directory and file of output

        Returns
        -------
        Bool
        """

        args = ["python", "../scripts/from_bed_to_BAM_for_chicago.py",
              bed, ncpus, output]

        logger.info("from_bed_to_BAM_for_chicago arguments:"+ " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        process.wait()

        if os.path.isfile(output + ".bam") is True:
            if os.path.getsize(output + ".bam") > 0:
                return True
        else:
            logger.fatal("from_bed_to_BAM_for_chicago\
                generates no output")
            return False

    def sort_output(self, output):
        """
        This function sort the output using samtools
        
        Parameters
        ----------
        output: str
            path to output directory and file
        """
        args = ["samtools", "sort",
                "-n", output+".bam",
                ">", output+"_sorted.bam"]

        logger.info("samtools args:"+ " ".join(args))

        process = subprocess.Popen(
            " ".join(args), shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        process.wait()
    
        if os.path.isfile(output+"_sorted.bam") is True:
            if os.path.getsize(output+"_sorted.bam") > 0:
                return True
        else:
            logger.fatal("samtools didnt generate sorted output")
            return False


    def run(self, input_files, input_metadata, output_files):
        """
        This function runs the wrapper_bed2chicago function
        and produce the output

        Parameters
        ----------
        input_files: dict
            deb
            ncpus
        input_metadata: dict
        output_files: dict
            output

        Returns
        -------
        results: bool
        output_metadata:dict
        """

        results = self.wrapper_bed2bam(
            input_files["bed"],
            input_files["ncpus"],
            output_files["output"])

        if results is True:
            sorted_results = self.sort_output(
                output_files["output"])

        output_metadata = {
            "BAM": Metadata(
                data_type="",
                file_type="bam",
                file_path=output_files["output"]+".bam",
                sources=input_metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                )
        }

        return sorted_results, output_metadata