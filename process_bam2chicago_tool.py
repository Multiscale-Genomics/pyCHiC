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

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from tool.bam2chicago_tool import bam2chicagoTool

################################################

class process_bam2chicago(Workflow):
    """
    This class creates .chinput files.
    input files for CHiCAGO
    """
    def __init__(self, configuration=None):
        """
        initiate the class

        Parameters
        ----------
        configuration: dict
            dictionary with parameters for different Tools, indicating
            how to run each of them
        """

        logger.info("Initialising process_bam2chicago")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that run the tools to create
        .chinput files

        Parameters:
        ----------
        input_files: dict
            BAM: str
                path to BAM files
            RMAP: str
                path to RMAP file
            BAITMAP: str
                path to BAITMAP file

        metadata: dict
            input metadata

        Returns:
        --------
        bool
        output_metadata: dict
            metadata for .chinput file
        """

        bam2chicago_caller = bam2chicagoTool(self.configuration)
        output_files_bam2chicago, output_metadata_bam2chicago = bam2chicago_caller.run(
            {
                "BAM" : input_files["BAM"],
                "RMAP" : input_files["RMAP"],
                "BAITMAP" : input_files["BAITMAP"]
            },
            {
                "BAM" : metadata["BAM"],
                "RMAP" : metadata["RMAP"],
                "BAITMAP" : metadata["BAITMAP"]
            },
            {
                "chrRMAP": output_files["chrRMAP"],
                "chrBAITMAP": output_files["chrBAITMAP"],
                "sample_name": output_files["sample_name"]
            }
        )

        out_path = output_files["sample_name"] + "/sampleout.chinput"

        if os.path.getsize(out_path) > 0:
            return output_files_bam2chicago, output_metadata_bam2chicago
        else:
            logger.fatal("process_bam2chicago failed to generate .chinput files")

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
    results = app.launch(process_bam2chicago,
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
        description="Pipeline to generate .chinput file")

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
        sys._run_from_cmdl = True # pylint: disable=protected-access

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
