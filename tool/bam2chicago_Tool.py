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
import shlex
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

    def check_chr_format(self, rmap_file, baitmap_file):
        """
        This function check that the chromosome format is right for
        bam3chicago.sh
        """
        rmapfile_new = rmap_file.split(".")[0]+"rfmt.rmap"
        baitmapfile_new = baitmap_file.split(".")[0]+"rfmt.baitmap"

        with open(rmap_file, "r") as file_in:
            for line in file_in:
                if line[0:3] == "chr":
                    rmapfile_new = rmap_file
                    break
                else:
                    line_hdl = line.rstrip().split("\t")
                    with open(rmapfile_new, "a") as file_out:
                        print("chr"+line, file=file_out, end="")
            logger.info("Checking rmap chromosome format..")


        with open(baitmap_file, "r") as file_in:
            for line in file_in:
                if line[0:3] == "chr":
                    baitmapfile_new = baitmap_file
                    break
                else:
                    line_hdl = line.rstrip().split("\t")
                    with open(baitmapfile_new, "a") as file_out:
                        print("chr"+line, file=file_out, end="")
            logger.info("Checking baitmap chromosome format..")


        return rmapfile_new, baitmapfile_new


    def bam2chicago(self, bamFiles, rmapFile, baitmapFile, sample_name):
        """
        Main function that preprocess the bam files into Chinput files. Part of
        the input files of CHiCAGO.
        It is a wrapper of bam2chicago.sh.

        Parameters
        ----------
        bamFiles : str,
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
        sample_name: str
            name of the output file. Bbam2chicago creates a folder with the
            name of this sample, and inside the folder there is a file with
            sample_name.chinput, that is the final output.



        Returns
        -------
        bool
        sample_name : str,
         name of the sample
        """

        args = ["../scripts/bam2chicago.sh",
                bamFiles,
                baitmapFile,
                rmapFile,
                sample_name]

        logger.info("bam2chicago CMD: " + " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            os.path.isfile(sample_name + ".chinput")
        except IOError:
            logger.fatal("bam2chicago failed to geenerate .chinput files")
            logger.fatal("bam2chicago stdout" + proc_out)
            logger.fatal("bam2chicago srerr"  + proc_err)
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Function that runs and pass the parameters to bam2chicago

        Parameters
        ----------
        input_files : dict
        bamfiles : str
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

        rfmat_rmap, rfmat_baitmap = self.check_chr_format(
            input_files["RMAP"],
            input_files["BAITMAP"])

        results = self.bam2chicago(
            input_files["BAM"],
            rfmat_rmap,
            rfmat_baitmap,
            output_files["sample_name"]
            )

        results = compss_wait_on(results)

        output_metadata = {
            "chinput": Metadata(
                data_type=input_metadata['BAM'].data_type,
                file_type="chinput",
                file_path=output_files["sample_name"],
                sources=[
                    input_metadata["BAM"].file_path,
                    input_metadata["RMAP"].file_path,
                    input_metadata["BAITMAP"].file_path
                ],
                taxon_id=input_metadata["BAM"].taxon_id,
                meta_data={
                    "assembly": input_metadata["BAM"].meta_data,
                    "tool": "bam2chicago"
                }
            )
        }

        return(results, output_metadata)


if __name__ == "__main__":

    path = os.path.join(os.path.dirname(__file__),"tests/data/")

    path = "/home/pablo/MuG/C-HiC/tests/data/"

    input_files = {
        "RMAP" : path + "test_makeRmap/test.rmap",
        "BAITMAP" : path +  "test_makeBaitmap/test.baitmap",
        "BAM" : path + "test_bed2bam/outbam_sorted.bam"
    }

    output_files = {
        "sample_name" :  path + "test_bam2chicago_Tool/sampleout"
    }

    input_metadata = {
        "RMAP" : Metadata(
            "data_chicago_input", ".rmap",
            path+"/h19_chr20and21_chr.rmap", None, {}, 9606),
        "BAITMAP" : Metadata(
            "data_chicago_input", ".baitmap",
            path+"/h19_chr20and21.baitmap_4col_chr.txt", None, {}, 9606),
        "BAM" : Metadata(
            "txt", "bamfile", path + "/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606)
    }

    bam2chicago_handle = bam2chicagoTool()
    bam2chicago_handle.run(input_files, input_metadata, output_files)

    out_path = output_files["sample_name"] + "/sampleout.chinput"

































