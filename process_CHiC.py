#!/usr/bin/env python

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

#Required for ReadTheDocs

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from CHiC.tool.truncater import Truncater
from CHiC.tool.rmap_tool import makeRmapFile
from CHiC.tool.makeBaitmap import makeBaitmapTool
from CHiC.tool.makeDesignFiles import makeDesignFilesTool
from CHiC.tool.fastq2bed import Fastq2bed
from CHiC.tool.bed2bam import bed2bam
from CHiC.tool.bam2chicago_tool import bam2chicagoTool
from CHiC.tool.run_chicago import ChicagoTool

################################################

class process_CHiC(Workflow):
    """
    This class output chromatin contacts from capture HiC from
    pair fastq reads, baits and RE
    """
    def __init__(self, configuration=None):
        """
        initiate the class

        Parameters
        ----------
        configuration: dict
            "RE_truncater" : str - format A^AGCT (Truncater),
            "RE" : {"HindIII" : 'A|AGCTT'},
        """

        logger.info("Initialising process_chicago_CHiC")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that run the tools to create
        .baitmap.

        Parameters:
        ----------
        input_files: dict
            fastq1: str
            fastq2: str
            genome_fa: str
                genome in fasta format

        input_metadata: dict
            input metadata

        output_files: dict

        Returns:
        --------


        #produce rmap file
        """

        try:
            rmap_caller = makeRmapFile(self.configuration)
            output_files_rmap, output_metadata_rmap = rmap_caller.run(
                {
                    "genome_fa" : input_files["genome_fa"]
                },
                {
                    "genome_fa" : metadata["genome_fa"]
                },
                {
                    "RMAP" : output_files["RMAP"],
                    "Rtree_file_dat" : output_files["Rtree_file_dat"],
                    "Rtree_file_idx" : output_files["Rtree_file_idx"]
                }
            )

            logger.info(".rmap file generated succesfully")

        except IOError:
            logger.fatal("rmap_tool failed to generate .rmap file")
            return False

        #produce Baitmap file
        try:
            baitmap_caller = makeBaitmapTool(self.configuration)
            output_files_baitmap, output_metadata_baitmap = baitmap_caller.run(
                {
                    "genome_idx" : input_files["genome_idx"],
                    "probes_fa": input_files["probes_fa"],
                    "Rtree_file_dat": input_files["Rtree_file_dat"],
                    "Rtree_file_idx": input_files["Rtree_file_idx"],
                    "genome_fa" : input_files["genome_fa"]
                },
                {
                    "genome_fa" : metadata["genome_fa"],
                    "probes_fa" : metadata["probes_fa"],
                    "Rtree_file_dat": metadata["Rtree_file_dat"],
                    "Rtree_file_idx": metadata["Rtree_file_idx"],
                    "genome_idx": metadata["genome_idx"]
                },
                {
                    "bait_sam" : output_files["bait_sam"],
                    "out_baitmap" : output_files["out_baitmap"],
                    "out_bam" : output_files["out_bam"]
                }
            )


            logger.info(".baitmap file generated succesfully")

        except IOError:
            logger.fatal("generate_CHiCAGO_baitmap failed to generate .baitmap file")
            return False


        #produce makeDesignFiles
        try:
            design_caller = makeDesignFilesTool(self.configuration)
            design_out, design_meta = design_caller.run(
                {
                    "RMAP" : input_files["RMAP"],
                    "BAITMAP": input_files["BAITMAP"]
                },
                {
                    "RMAP" : metadata["RMAP"],
                    "BAITMAP" : metadata["BAITMAP"]
                },
                {
                    ".nbpb" : output_files[".nbpb"],
                    ".npb"  : output_files[".npb"],
                    ".poe" : output_files[".poe"]
                }
            )

            logger.info("Design files succesfully generated")

        except IOError:
            logger.fatal("process_makeDesign failed to" +
                         "generate design files")
            return False

        try:
            truncater_caller = Truncater(self.configuration)
            output_files_truncater, output_metadata_truncater = truncater_caller.run(
                {
                    "fastq1": input_files["fastq1"],
                    "fastq2": input_files["fastq2"]
                },
                {
                    "fastq1": metadata["fastq1"],
                    "fastq2": metadata["fastq2"]
                },
                {
                    "fastq1_trunc": output_files["fastq1_trunc"],
                    "fastq2_trunc" : output_files["fastq2_trunc"],
                    "hicup_summary" : output_files["hicup_summary"],
                    "barchat_fastq1" : output_files["barchat_fastq1"],
                    "barchat_fastq2" : output_files["barchat_fastq2"]
                }
            )

        except IOError:
            return False

        output_files = {}
        output_metadata = {}

        output_files.update(output_files_rmap)
        output_files.update(output_files_baitmap)
        output_files.update(design_out)
        output_files.update(output_files_truncater)

        output_metadata.update(output_metadata_rmap)
        output_metadata.update(output_metadata_baitmap)
        output_metadata.update(design_meta)
        output_metadata.update(output_files_truncater)

        return output_files, output_metadata

#############################################################

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function

    This function launch the app using the configuration written
    in two json files: config_process_rmapBaitmap.json  and
    input_process_rmapBaitmap.json
    """
    #Instantiate and lauch the app
    print("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    results = app.launch(process_CHiC,
                         config,
                         in_metadata,
                         out_metadata)
    #2. The App has finished
    print("2. Execution finished: see " + out_metadata)
    print(results)

    return results

#########################################################################


if __name__ == "__main__":

    #set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Pipeline to generate .baitmap file")

    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument(
        "--in_metadata", help="Location of metadata file")
    PARSER.add_argument(
        "--out_metadata", help="Location of output metadata file")
    PARSER.add_argument(
        "--local", action="store_const", const=True, default=False)

    #Get matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata
    LOCAL = ARGS.local

    if LOCAL:
        import sys
        sys._run_from_cmdl = True # pylint: disable=protected-access

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
