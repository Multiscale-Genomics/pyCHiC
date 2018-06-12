
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
    This class contain functions to convert a bed file bam_out from
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
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def wrapper_bed2bam(self, bed, bam_out):
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
        bam_out: str
            path to the bam_out directory and file of bam_out

        Returns
        -------
        Bool
        """

        args = ["python", "../scripts/from_bed_to_bam.py",
                bed, "2", bam_out]

        logger.info("from_bed_to_BAM_for_chicago arguments:"+ " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        if os.path.isfile(bam_out + ".bam") is True:
            if os.path.getsize(bam_out + ".bam") > 0:
                return True
        else:
            logger.fatal("from_bed_to_BAM_for_chicago\
                generates no bam_out")
            return False

    def sort_bam_out(self, bam_out):
        """
        This function sort the bam_out using samtools

        Parameters
        ----------
        bam_out: str
            path to bam_out directory and file
        """
        args = ["samtools", "sort",
                "-n", bam_out+".bam",
                ">", bam_out+"_sorted.bam"]

        logger.info("samtools args:"+ " ".join(args))

        process = subprocess.Popen(
            " ".join(args), shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        process.wait()

        if os.path.isfile(bam_out+"_sorted.bam") is True:
            if os.path.getsize(bam_out+"_sorted.bam") > 0:
                return True
        else:
            logger.fatal("samtools didnt generate sorted bam_out")
            return False

    def run(self, input_files, input_metadata, output_files):
        """
        This function runs the wrapper_bed2chicago function
        and produce the bam_out

        Parameters
        ----------
        input_files: dict
            deb
            ncpus
        input_metadata: dict
        output_files: dict
            bam_out

        Returns
        -------
        results: bool
        output_metadata:dict
        """

        results = self.wrapper_bed2bam(
            input_files["bed"],
            output_files["bam_out"])

        if results is True:
            sorted_results = self.sort_bam_out(
                output_files["bam_out"])

        output_metadata = {
            "BAM": Metadata(
                data_type="",
                file_type="bam",
                file_path=output_files["bam_out"]+".bam",
                sources=input_metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                )
        }

        return sorted_results, output_metadata
