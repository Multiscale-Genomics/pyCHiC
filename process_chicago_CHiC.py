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
import os

#Required for ReadTheDocs
from functools import wraps

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from tool.truncater import Truncater
from tool.makeRmap_tool import makeRmapFile
from tool.makeBaitmap import makeBaitmapTool
from tool.makeDesignFiles_Tool import makeDesignFilesTool
from tool.fastq2bed import Fastq2bed
from tool.bed2bam import bed2bam
from tool.bam2chicago_Tool import bam2chicagoTool
from tool.runChicago import ChicagoTool

################################################

class process_chicago_CHiC(Workflow):
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

    def run(self, input_files, input_metadata, output_files):
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

        """
        #truncate reads using HiCUP truncater
        truncater_caller = Truncater(self.configuration)
        output_files_Truncater, output_metadata_Truncater = truncater_caller.run(
            {
                "fastq1" : input_files["fastq1"],
                "fastq2": input_files["fastq2"],
            },
            {
                "fastq1" : input_metadata["fastq1"],
                "fastq2" : input_metadata["fastq2"],
            },
            {
                "out_dir" : output_files["out_dir"],
            }
        )

        fastq1_name = input_files["fastq1"].split("/")[-1].split(".")[-2]+".trunc.fastq"
        fastq2_name = input_files["fastq2"].split("/")[-1].split(".")[-2]+".trunc.fastq"

        if os.path.getsize(output_files["out_dir"]+fastq1_name) > 0 and \
           os.path.getsize(output_files["out_dir"]+fastq2_name) > 0:
            pass

        else:
            logger.fatal("Truncater failed to truncate  fastq reads")
            return False


        #make rmap and rtree files

        rmap_caller = makeRmapFile(self.configuration)
        output_files_rmap, output_meta_rmap = rmap_caller.run(
            {
                "genome_fa": input_files["genome_fa"],
            }, {
                "genome_fa": input_metadata["genome_fa"]
            }, {
                "out_dir_makeRmap": output_files["out_dir"],
                "out_prefix_makeRmap": output_files["out_prefix_makeRmap"],
                "Rtree_files": output_files["Rtree_files"]
            }
        )

        out = "".join(
            [
                f for f in os.listdir(output_files["out_dir"])
                if f.startswith("Digest_") and f.endswith(".map")
            ]
        )

        if os.path.getsize(output_files["out_dir"] + out) > 0 and \
           os.path.getsize(output_files["Rtree_files"] + ".dat") and \
           os.path.getsize(output_files["Rtree_files"] + ".idx"):
            pass
        else:
            logger.fatal("makeRmapFile failed to generate rmap file")
            return False

        #make baitmap file

        baitmap_called = makeBaitmapTool(self.configuration)
        output_files_baitmap, output_meta_baitmap = baitmap_called.run(
            {

            }
            )




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
    results = app.launch(generate_CHiCAGO_baitmap,
                         config,
                         in_metadata,
                         out_metadata)
    #2. The App has finished
    print("2. Execution finished: see " + out_metadata)
    print(results)

    return results

#########################################################################


if __name__ == "__name__":

    #set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Pipeline to generate .baitmap file")

    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument(
        "--in_metadata", help="Location of metadata file")
    PARSER.add_argument(
        "--out_metadata", help="Location of output metadata file")
    PARSER.add_argument(
        "--local", action="store_const", cont=True, default=False)

    #Get matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata
    LOCAL = ARGS.local

    if LOCAL:
        import sys
        sys._run_from_cmdl = True

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
