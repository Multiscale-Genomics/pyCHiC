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

    @staticmethod
    def check_chr_format(rmap_file, baitmap_file, chrRMAP, chrBAITMAP):
        """
        This function check that the chromosome format is right for
        bam2chicago.sh
        """
        rmapfile_new = rmap_file.split(".")[0]+"rfmt.rmap"
        baitmapfile_new = baitmap_file.split(".")[0]+"rfmt.baitmap"
        rformat_rmap = True
        rformat_baitmap = True

        logger.info("checking .rmap chr format")
        try:
            rmapfile_new = pd.read_csv(rmap_file, header=None, sep="\t")
        except IOError:
            logger.fatal("bam2chicago_tool failed on checking the rmap format 1")
	    print(rmap_file)
        try:
            chr_rmap = str(rmapfile_new.iloc[0, 0])
        except:
            logger.fatal("bam2chicago_tool failed on checking the rmap format 2")
        try:
            if str(chr_rmap[0:3]) == "chr":
                logger.info("rmap file chromosome on the right format")
                rformat_rmap = False
            else:
                rmapfile_new.iloc[:, 0] = rmapfile_new.iloc[:, 0].apply(lambda x: "chr"+str(x))
        except:
            logger.fatal("bam2chicago_tool failed on checking the rmap format 3")

        logger.info("checking .baitmap chr format")

        try:
            baitmapfile_new = pd.read_table(baitmap_file, header=None, sep="\t")
            chr_baitmap = str(baitmapfile_new.iloc[0, 0])
            if chr_baitmap[0:3] == "chr":
                logger.info("baitmap file chromosome on the right format")
                rformat_baitmap = False

            else:
                baitmapfile_new.iloc[:, 0] = baitmapfile_new.iloc[:, 0].apply(lambda x: "chr"+str(x))

        except IOError:
            logger.fatal("bam2chicago_tool failed on checking the baitmap format")


        if rformat_rmap:
            rmapfile_new.to_csv(chrRMAP, sep="\t", header=None, index=False)
        else:
            chrRMAP = rmap_file

        if rformat_baitmap:
            baitmapfile_new.to_csv(chrBAITMAP, sep="\t", header=None, index=False)
        else:
            chrBAITMAP = baitmap_file

        return chrRMAP, chrBAITMAP

    @staticmethod
    def bam2chicago(bamFile, rmapFile, baitmapFile, chinput):
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
        no_tar_out = "".join(chinput.split(".")[0])

        try:
            bam2chicago_script = os.path.join(os.path.dirname(__file__), "scripts/bam2chicago.sh")

            args = [bam2chicago_script,
                    bamFile,
                    baitmapFile,
                    rmapFile,
                    no_tar_out]

            logger.info("bam2chicago CMD: " + " ".join(args))

            process = subprocess.Popen(
                ' '.join(args),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

            process.wait()

            try:
                common.tar_folder(
                    no_tar_out,
                    chinput,
                    os.path.split(no_tar_out)[1]
                    )

            except IOError:
                logger.fatal("could not tar folder")
                print(no_tar_out, chinput,os.path.split(no_tar_out)[1])

            return True
        except IOError:
            logger.fatal("bam2chicago failed to generate chicago output files =(")
            return False


    def run(self, input_files, metadata, output_files):
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
        if os.path.isdir(os.path.split(output_files["chinput"])[0]) is False:
            logger.info("creating output directory")
            os.mkdir(os.path.split(output_files["chinput"])[0])

        folder_name = os.path.split(input_files["hicup_outdir_tar"])[0] + "/"+\
                    "".join(os.path.split(input_files["hicup_outdir_tar"])[1].split(".")[:-1])

        tar = tarfile.open(input_files["hicup_outdir_tar"])
        tar.extractall(path="".join(os.path.split(input_files["hicup_outdir_tar"])[0]))
        tar.close()

        bam_file = "".join([file_hdl for file_hdl in os.listdir(folder_name)
                            if file_hdl.endswith(".bam")])

        path_bam =  folder_name + "/" + bam_file

        if self.configuration["aligner"] == "tadbit":
            logger.info("cheking chr format from rmap and baitmap")
            rfmat_rmap, rfmat_baitmap = self.check_chr_format(
                input_files["RMAP"],
                input_files["BAITMAP"],
                output_files["chrRMAP"],
                output_files["chrBAITMAP"])

        else:
            rfmat_rmap = input_files["RMAP"]
            rfmat_baitmap = input_files["BAITMAP"]

        results = self.bam2chicago(
            path_bam,
            rfmat_rmap,
            rfmat_baitmap,
            output_files["chinput"]
            )

        results = compss_wait_on(results)

        output_metadata = {
            "chinput" : Metadata(
                data_type="CHiC_data",
                file_type="tar",
                file_path=output_files["chinput"],
                sources=[
                    metadata["RMAP"].file_path,
                    metadata["BAITMAP"].file_path,
                    metadata["hicup_outdir_tar"].file_path
                    ],
                taxon_id=metadata["hicup_outdir_tar"].taxon_id,
                meta_data={"tool": "bam2chicago_tool"}
            ),
            "chrRMAP" : Metadata(
                data_type=metadata['RMAP'].data_type,
                file_type=".rmap",
                file_path=rfmat_rmap,
                sources=[
                    metadata["RMAP"].file_path,
                    ],
                taxon_id=metadata["RMAP"].taxon_id,
                meta_data={"tool": "bam2chicago_tool"}
            ),
            "chrBAITMAP" : Metadata(
                data_type=metadata['BAITMAP'].data_type,
                file_type=".baitmap",
                file_path=rfmat_baitmap,
                sources=[
                    metadata["BAITMAP"].file_path
                    ],
                taxon_id=metadata["BAITMAP"].taxon_id,
                meta_data={"tool": "bam2chicago_tool"}
            )
        }

        return output_files, output_metadata
