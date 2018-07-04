
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

        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, fastq1=FILE_IN, fastq2=FILE_IN,
          fastq1_trunc=FILE_OUT, fastq2_trunc=FILE_OUT,
          parameters=IN)
    def truncate_reads(self, fastq1, fastq2, fastq1_trunc, fastq2_trunc,
                       parameters):
        """
        Function to truncate the reads with hicup_truncater

        Parameters
        ----------
        fastq1: str
            path to fastq file read 1
        fastq2: str
            path to feastq file read 2
        fastq1_trunc:
            path to the fastq1_trunc file
            IMP: This file shoudl be called as the input file. with
                Name.<trunc>.fastq
        fastq2_trunc:
            path to the fastq1_trunc file
            IMP: This file shoudl be called as the input file. with
                Name.<trunc>.fastq
        parameters: list
            paramters selected using get_params()

        Return
        ------
        bool
        """
        temp_fastq1 = "".join(fastq1.split("/")[-1])
        temp_fastq2 = "".join(fastq2.split("/")[-1])
        temp_fastq1_trunc = "".join(fastq1_trunc.split("/")[-1])
        temp_fastq2_trunc = "".join(fastq2_trunc.split("/")[-1])

        copy(fastq1, temp_fastq1)
        copy(fastq1, temp_fastq2)

        args = ["hicup_truncater", temp_fastq1, temp_fastq2]

        args += parameters

        logger.info("hicup_truncater command: "+ " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        copy(temp_fastq1_trunc, fastq1_trunc)
        copy(temp_fastq2_trunc, fastq2_trunc)

        os.remove(temp_fastq1)
        os.remove(temp_fastq2)
        os.remove(temp_fastq1_trunc)
        os.remove(temp_fastq2_trunc)

        if os.path.isfile(fastq1_trunc) is True:
            if os.path.getsize(fastq1_trunc) > 0:
                if os.path.isfile(fastq2_trunc) is True:
                    if os.path.getsize(fastq2_trunc) > 0:
                        return True

        return False

    @staticmethod
    def get_params(configuration):
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
            "outdir": ["--outdir", True]
        }

        params = []

        for arg in configuration:
            if arg in parameters:
                if parameters[arg][1] is True:
                    if parameters[arg][0] == "--outdir":
                        name = "."
                        params += [parameters[arg][0], name]
                    else:
                        params = [parameters[arg][0], configuration[arg]]
                else:
                    params += [parameters[arg][0]]

        return params

    def run(self, input_files, metadata, output_files):
        """
        This is the function that runs all the functions

        Parameters
        ----------
        input_files: dict
            fastq1
            fastq2
        metadata: dict
        output_files: dict
            out_dir: str
                directory to write the output

        Return
        ------
        output_files: dict
            fastq1_trunc
            fastq2_trunc
        output_metadata: dict
        """

        param_truncater = self.get_params(self.configuration)

        out_fastq1_trunc = "/".join(output_files["fastq1_trunc"].split("/")[:-1])
        out_fastq2_trunc = "/".join(output_files["fastq2_trunc"].split("/")[:-1])

        if out_fastq1_trunc != out_fastq2_trunc:
            logger.fatal("Output directory for the truncated reads shoudl be "+
                         "the same", out_fastq1_trunc, out_fastq2_trunc)

        if out_fastq1_trunc != self.configuration["outdir"]:
            logger.fatal("outdir parameter should be the same directory as the fastq1_trunc"+
                         " and fastq2_trunc: "+ out_fastq1_trunc + " " +
                         self.configuration["outdir"])

        out_dir = self.configuration["outdir"]

        logger.info("truncater parameters: "+ param_truncater)

        results = self.truncate_reads(
            input_files["fastq1"],
            input_files["fastq2"],
            output_files["fastq1_trunc"],
            output_files["fastq2_trunc"],
            param_truncater)

        results = compss_wait_on(results)

        output_metadata = {
            "fastq1_trunc": Metadata(
                data_type="FASTQ",
                file_type="FASTQ",
                file_path=output_files["fastq1_trunc"],
                sources=[metadata["fastq1"].file_path, metadata["fastq2"].file_path],
                taxon_id=9606,
                meta_data=""
            ),
            "fastq2_trunc": Metadata(
                data_type="text",
                file_type="tsv",
                file_path=output_files["fastq1_trunc"],
                sources=[metadata["fastq1"].file_path, metadata["fastq2"].file_path],
                taxon_id=9606,
                meta_data=""
            )
        }

        return output_files, output_metadata
