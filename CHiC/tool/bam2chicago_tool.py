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

import os
import subprocess
import sys
from shutil import move
from shutil import rmtree
import tarfile
import pandas as pd
from tool.common import common
from utils import logger
import re

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

#######################################################

class bam2chicagoTool(Tool):
    """
    Tool for preprocess the input files
    """

    def __init__(self, configuration=None):
        """
        initialising the function

        Parameters:
        -----------
        configuration: dict
         dictionary containing all the arguments and parameters
         to run the tool
        """

        print("bam2chicago initialising")
        Tool.__init__(self)


        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, bamFile=FILE_IN, rmapFile=FILE_IN, baitmapFile=FILE_IN,
          chinput=FILE_OUT)
    def bam2chicago(self, bamFile, rmapFile, baitmapFile, chinput):
        """
        Main function that preprocess the bam files into Chinput files. Part of
        the input files of CHiCAGO.
        It is a wrapper of bam2chicago.sh.

        Parameters
        ----------
        bamFile : str,
            path to paired-end file produced by a HiC aligner; Chicago has
            only been tested with data produced by HiCUP
            (http://www.bioinformatics.babraham.ac.uk/projects/hicup/).
            However, it should theoretically be possible to use other HiC
            aligners for this purpose.
        rmapFile : str,
            A tab-separated file of the format
            <chr> <start> <end> <numeric ID>,
            describing the restriction digest (or "virtual digest"
            if pooled fragments are used). These numeric IDs are referred to as
            "otherEndID" in Chicago. All fragments mapping outside of the digest
            coordinates will be disregarded by both these scripts and Chicago.
        baitMapFile: str,
            Tab-separated file of the format
            <chr> <start> <end> <numeric ID> <annotation>,
            listing the coordinates of the baited/captured
            restriction fragments (should be a subset of the fragments
            listed in rmapfile), their numeric IDs (should match those listed
            in rmapfile for the corresponding fragments) and their annotations
            (such as, for example, the names of baited promoters). The numeric
            IDs are referred to as "baitID" in Chicago.
        chinput: str
            name of the output file. Bbam2chicago creates a folder with the
            name of this sample, and inside the folder there is a file with
            chinput.chinput, that is the final output.



        Returns
        -------
        bool
        chinput : str,
         name of the sample
        """
        out_folder = "".join(chinput.split(".")[0])



        try:
            bam2chicago_script = os.path.join(os.path.dirname(__file__), "scripts/bam2chicago.sh")

            args = [bam2chicago_script,
                    bamFile,
                    baitmapFile,
                    rmapFile,
                    out_folder]

            logger.info("bam2chicago CMD: " + " ".join(args))

            process = subprocess.Popen(
                ' '.join(args),
                shell=True
                )

            process.wait()

            path_out = out_folder+"/"+os.path.split(out_folder)[1]+".chinput"
            move(path_out, chinput)
            print("path_out "+path_out, "chinput "+chinput)

            return True
        except IOError:
            logger.fatal("bam2chicago failed to generate chicago output files =(")
            return False

    @task(returns=bool,
          sorted_bam=FILE_OUT,
          hicup_outdir_tar=FILE_IN,
          bam_name=IN
         )
    def sort_chicago(self, sorted_bam, hicup_outdir_tar, bam_name):
        """
        This function sort bamfile by name of the reads as bam2chicago requires

        Parameters
        ----------
        bamfile: str
        sorted_bam: str

        Returns
        -------
        sorted_bam: str
        """

        bamfile = self.untar_hicup_out(hicup_outdir_tar, bam_name)

        logger.info("Sorting bamfile")

        args = ["samtools", "sort", "-n", bamfile, "-o", sorted_bam]

        logger.info("samtools CMD: " + " ".join(args))

        process = subprocess.Popen(
            ' '.join(args),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        process.wait()

        return True

    @staticmethod
    def untar_hicup_out(hicup_outdir_tar, bam_name):
        """
        Untar hicup output filder

        Parameters
        ----------
        hicup_outdir_tar: str
            path to hicup output folder
        path_bam: str
            path to bam file

        Returns
        -------
        bool
        """
        logger.info("UNTAR: "+hicup_outdir_tar)
        logger.info("UNTAR PATH: "+os.path.split(hicup_outdir_tar)[0])
        try:
            tar = tarfile.open(hicup_outdir_tar)
            tar.extractall(path="".join(os.path.split(hicup_outdir_tar)[0]))
            tar.close()
            logger.info("bam2chicago succesfully untar hicup_outdir_tar")

            bamfile = os.path.split(hicup_outdir_tar)[0] + "/"+\
                    "".join(os.path.split(hicup_outdir_tar)[1].split(".")[:-1])+\
                    "/"+bam_name

        except IOError:
            logger.info("bam2chicago could not extract hicup_outdir_tar")

        logger.info("Bamfile: "+bamfile)
        return bamfile

    def run(self, input_files, input_metadata, output_files):
        """
        Function that runs and pass the parameters to bam2chicago

        Parameters
        ----------
        input_files : dict
        hicup_outdir_tar : str
        rmapFile : str
        baitmapFile : str

        metadata : dict

        Returns
        -------
        output_files : list
        List of locations for the output files.
        output_metadata : list
        List of matching metadata dict objects
        """
        RMAP = "tests/data/test_run_chicago/test.rmap"
        BAITMAP = "tests/data/test_run_chicago/test.baitmap"

        #hicup_outdir_tar = "tests/data/test_hicup/output.tar"
        output_files["hicup_outdir_tar"] = self.configuration["execution"]+"/"+\
                                           os.path.split(output_files["hicup_outdir_tar"])[1]


        output_files["chinput"] = self.configuration["execution"]+"/"+\
                                    os.path.split(output_files["chinput"])[1]

        #find the name of the future bamfile
        fastq1 = str(os.path.split(input_files["fastq1"])[1])
        fastq2 = str(os.path.split(input_files["fastq2"])[1])

        fastq1_1 = [pos for pos, char in enumerate(fastq1) if char == "1"]
        fastq2_2 = [pos for pos, char in enumerate(fastq2) if char == "2"]

        sample_indicator_index = list(set(fastq1_1).intersection(fastq2_2))
        sample_indicator_index = int("".join([str(n) for n in sample_indicator_index]))
        fastq1 = "".join(fastq1.split(".")[0])

        bam_name = fastq1[:sample_indicator_index]+ "1_2" + \
                   fastq1[sample_indicator_index+1:]+ ".hicup.bam"


        folder_name = os.path.split(output_files["hicup_outdir_tar"])[0] + "/"+\
                    "".join(os.path.split(output_files["hicup_outdir_tar"])[1].split(".")[:-1])

        #path_bam = folder_name + "/" + bam_file

        #bamfile = self.untar_hicup_out(output_files["hicup_outdir_tar"], bam_name)

        sorted_bam = self.configuration["execution"] + "/" + "sorted_bam"

        self.sort_chicago(sorted_bam,
                          output_files["hicup_outdir_tar"],
                          bam_name
                         )

        results = self.bam2chicago(
            sorted_bam,
            RMAP,
            BAITMAP,
            output_files["chinput"]
            )

        #rmtree(folder_name)

        output_metadata = {
            "chinput" : Metadata(
                data_type="data_chic",
                file_type="TXT",
                file_path=output_files["chinput"],
                sources=[
                    RMAP,
                    BAITMAP,
                    sorted_bam
                    ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={"tool": "process_CHiC",
                           "tool_description" : "bam2chicago_tool"}
            )
        }

        return output_files, output_metadata
