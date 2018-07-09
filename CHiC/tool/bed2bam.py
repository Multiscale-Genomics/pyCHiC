
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
from shutil import copy

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

    @task(returns=bool, bed=FILE_IN, bam_out=FILE_OUT)
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
        script = os.path.join(os.path.dirname(__file__), "scripts/from_bed_to_bam.py")

        print(os.getcwd())

        args = ["python", script,
                bed+"_temp", self.configuration["ncpus"], "outbam_temp"]

        logger.info("from_bed_to_BAM_for_chicago arguments:"+ " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        try:
            with open("outbam_temp", "w") as f_out:
                process = subprocess.Popen(
                    ' '.join(args),
                    shell=True,
                    stdout=f_out, stderr=f_out
                    )
                process.wait()
            return True

        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, args))
            return False

        try:
            with open("outbam_temp", "r") as f_in:
                with open(bam_out, "w") as f_out:
                    f_out.write(f_in.read())
            return True

        except IOError:
            return False


    @task(returns=bool, bam_out=FILE_IN,
          bam_out_sorted=FILE_OUT)
    def sort_bam_out(self, bam_out, bam_out_sorted):
        """
        This function sort the bam_out using samtools

        Parameters
        ----------
        bam_out: str
            path to bam_out directory and file
        """
        args = ["samtools", "sort",
                "-n", bam_out]

        logger.info("samtools args:"+ " ".join(args))

        try:
            with open(bam_out_sorted, "w") as f_out:
                process = subprocess.Popen(
                    ' '.join(args),
                    shell=True,
                    stdout=f_out, stderr=f_out
                    )
                process.wait()
            return True

        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, args))
            return False

    def run(self, input_files, metadata, output_files):
        """
        This function runs the wrapper_bed2chicago function
        and produce the bam_out

        Parameters
        ----------
        input_files: dict
            deb
            ncpus
        metadata: dict
        output_files: dict
            bam_out

        Returns
        -------
        results: bool
        output_metadata:dict
        """
        output_dir = "/".join(output_files["bam_out"].split("/")[:-1])

        if os.path.isdir(output_dir) is False:
            logger.info("creating output directory")
            os.mkdir(output_dir)

        results = self.wrapper_bed2bam(
            input_files["bed"],
            output_files["bam_out"])

        results = compss_wait_on(results)

        if results is True:
            sorted_results = self.sort_bam_out(
                output_files["bam_out"], output_files["bam_out_sorted"])

        output_metadata = {
            "bam_out": Metadata(
                data_type="TXT",
                file_type="bam",
                file_path=output_files["bam_out"],
                sources=metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                ),
            "bam_out_sorted": Metadata(
                data_type="TXT",
                file_type="bam",
                file_path=output_files["bam_out_sorted"],
                sources=metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                ),
        }

        return output_files, output_metadata
