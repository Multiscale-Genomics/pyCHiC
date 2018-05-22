
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


class Fastq2bed(Tool):
    """
    This class has the functions to convert fastq
    files to bam files.
    Hicup_truncated, gem(binaries) and tadbit is neccesary to be
    installed, and added to $PATH.

    """
    def __init__(self, configuration=None):
        """
        initialising the clas

        Parameters
        ----------
        configuration: dict
            contains parameters to run the functino
        """
        logger.info("initialising Fastq2bed")

    def tadbit_map(self, fastq1, fastq2, genindex, RE, wd, chromosome):
        """
        This function map the Capture fastq reads to the reference genome using gem

        Parameters
        ----------
        fastq1: str
            path to fastq
        fastq2: str
            path to fastq
        genindex: str
            path to the ref. genome indexed with gem.
        RE: str
            restriction enzyme used to digest the genome.
            the form of the RE is the name with case
            sensitive. exam. HindIII
        wd: str
            working diretory for the output
        chr: str
            chromosome number to run just one chromosome.
            chr number

        Returns
        -------
        bool
        """
        logger.info(chromosome)
        if chromosome is "":
            args1 = ["tadbit", "map",
                     "--fastq", fastq1,
                     "--index", genindex,
                     "--read", "1",
                     "--renz", RE,
                     "-w", wd]

            args2 = ["tadbit", "map",
                     "--fastq", fastq2,
                     "--index", genindex,
                     "--read", "2",
                     "--renz", RE,
                     "-w", wd]
        else:
            args1 = ["tadbit", "map",
                     "--fastq", fastq1,
                     "--index", genindex,
                     "--read", "1",
                     "--renz", RE,
                     "-w", wd,
                    "--chr_name", chromosome]

            args2 = ["tadbit", "map",
                     "--fastq", fastq2,
                     "--index", genindex,
                     "--read", "2",
                     "--renz", RE,
                     "-w", wd,
                     "--chr_name", chromosome]


        logger.info("tadbit map1 arguments:" + " ".join(args1))

        process1 = subprocess.Popen(" ".join(args1), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process1.wait()

        try:
            mapped_files = os.listdir(wd+"/01_mapped_r1/")
            for indv_file in mapped_files:
                if "full" in indv_file:
                    if os.path.getsize(wd+"/01_mapped_r1/"+indv_file) > 0:
                        pass
        except IOError:
            logger.fatal("tadbit map1 failed to generate" \
                          "mapped files")
            return False

        logger.info("tadbit map2 arguments:"+ " ".join(args2))

        process2 = subprocess.Popen(" ".join(args2), shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        process2.wait()

        try:
            mapped_files = os.listdir(wd+"/01_mapped_r2/")
            for indv_file in mapped_files:
                if os.path.getsize(wd+"/01_mapped_r2/"+indv_file) > 0:
                    pass
        except IOError:
            logger.fatal("tadbit map2 failed to generate" \
                          "mapped files")
            return False

        return True

    def tadbit_parse(self, genome_fasta, wd):
        """
        This function parse the output files from mapping

        Parameters
        ----------
        genome: str
            path to the reference genome in fasta
        wd: str
            working directory path

        Returns
        -------
        bool
        """

        args = ["tadbit", "parse",
                "--genome", genome_fasta,
                "-w", wd]

        logger.info("arguments for tadbit parse:" + " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        try:
            mapped_files = os.listdir(wd+"/02_parsed_reads/")
            for indv_file in mapped_files:
                if os.path.getsize(wd+"/02_parsed_reads/"+indv_file) > 0:
                    pass
        except IOError:
            logger.fatal("tadbit parse failed to generate" \
                          "mapped files")
            return False

        return True

    def tadbit_filter(self, wd):
        """
        This function filter the reads

        Parameters
        ----------
        wd: str
            path to the workign directory

        Returns
        -------
        bool
        """

        args = ["tadbit", "filter",
                "-w", wd]

        logger.info("arguments for tadbit filter:" + " ".join(args))

        process = subprocess.Popen(" ".join(args), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()

        try:
            mapped_files = os.listdir(wd+"/03_filtered_reads/")
            for indv_file in mapped_files:
                if os.path.getsize(wd+"/03_filtered_reads/"+indv_file) > 0:
                    pass
        except IOError:
            logger.fatal("tadbit parse failed to generate" \
                          "mapped files")
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        This function run all the functions and generate the output files

        Parameters
        ----------
        input_files: dict
            fastq1,
            fastq2,
            RE,
            chromosome,
            genindex
            genome_fasta
        input_metadata: dict
        output_files: dict
            wd

        Returns
        -------
        bool
        """
        output_metadata = {
            "deb": Metadata(
                data_type="text",
                file_type="tsv",
                file_path=output_files["wd"]+"/03_filtered_reads",
                sources=[input_metadata["fastq1"].file_path,
                         input_metadata["fastq2"].file_path],
                taxon_id=9606,
                meta_data=""
            )
        }

        results_map = self.tadbit_map(input_files["fastq1"],
                                      input_files["fastq2"],
                                      input_files["genindex"],
                                      input_files["RE"],
                                      output_files["wd"],
                                      input_files["chromosome"])

        if results_map is True:
            results_parse = self.tadbit_parse(input_files["genome_fasta"],
                                              output_files["wd"])
            if results_parse is True:
                results_filter = self.tadbit_filter(output_files["wd"])

                return results_filter, output_metadata

        return False
