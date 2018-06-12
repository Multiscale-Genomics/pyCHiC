
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


class Truncater(Tool):
    """
    This class has the functions to truncate fastq reads with
    Hicup_truncated. Hicup must be installed and all files in the
    PATH
    """
    def __init__(self, configuration=None):
        """
        initialising the class

        Parameters
        ----------
        configuration: dict
            contains parameters to run the functino
        """
        logger.info("Initialising truncater")
        self.configuration.update(configuration)

    def truncate_reads(self, fastq1, fastq2, parameters, out_dir):
        """
        Function to truncate the reads with hicup_truncater

        Parameters
        ----------
        fastq1: str
            path to fastq file read 1
        fastq2: str
            path to feastq file read 2
        parameters: list
            paramters selected using get_params()

        Return
        ------
        bool
        """
        args = ["hicup_truncater", fastq1, fastq2, "--outdir", out_dir,
                "--quiet"]

        args += parameters

        logger.info("hicup_truncater parameters: "+ " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        name_trunc1 = ""
        name_trunc2 = ""
        if "/" in fastq1:
            name_trunc1 = fastq1.split("/")[-1]
            name_trunc1 = name_trunc1.split(".")[-2]+".trunc.fastq"
        else:
            name_trunc1 = name_trunc1.split(".")[-2]+".trunc.fastq"

        if "/" in fastq2:
            name_trunc2 = fastq2.split("/")[-1]
            name_trunc2 = name_trunc2.split(".")[-2]+".trunc.fastq"
        else:
            name_trunc2 = name_trunc2.split(".")[-2]+".trunc.fastq"

        if os.path.isfile(out_dir+name_trunc1) is True:
            pass
        else:
            return False

        if os.path.getsize(out_dir+name_trunc1) > 0:
            return True

        return False

    def get_params(self,configuration):
        """
        this function take the parameters that
        have been selected tu run the truncater

        Paramaters
        ----------
         --outdir: str,
            output directory
         --quiet: str, flag
            supress all the progress report
        --re1: str
            restriction enzyme to recognise. A^AGCT
        --sequences: str
            instead of RE, recognise a DNA sequence
        --threads: int
            number of threads to use
        --zip: str, flag
            compress output

        Returns
        -------
        arguments: list
           list with arguments to pass to truncate_read function
        """

        parameters = {
            "quiet_progress": ["--quiet", False],
            "RE_truncater": ["--re1", True],
            "sequence_junction": ["--sequences", True],
            "threads" : ["--threads", True],
            "zip": ["--zip", False],
        }

        params = []

        for arg in configuration:
            if arg in parameters:
                if parameters[arg][1] is True:
                    params += [parameters[arg][0], configuration[arg]]
                else:
                    params += [parameters[arg][0]]

        return params

    def run(self, input_files, input_metadata, output_files):
        """
        This is the function that runs all the functions

        Parameters
        ----------
        input_files: dict
            fastq1
            fastq2
        input_metadata: dict
        output_files: dict
            out_dir: str
                directory to write the output

        Return
        ------
        output_files: dict
            fastq_trunc1
            fastq_trunc2
        output_metadata: dict
        """

        param_truncater = self.get_params(self.configuration)

        results = self.truncate_reads(
            input_files["fastq1"],
            input_files["fastq2"],
            param_truncater,
            output_files["out_dir"])

        output_metadata = {
            "tsv": Metadata(
                data_type="text",
                file_type="tsv",
                file_path=output_files["out_dir"],
                sources=[input_metadata["fastq1"].file_path, input_metadata["fastq2"].file_path],
                taxon_id=9606,
                meta_data=""
            )
        }

        return results, output_metadata