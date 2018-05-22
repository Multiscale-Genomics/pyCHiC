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

from tool.makeRmap_tool import makeRmapFile

################################################

class generate_CHiCAGO_rmap(Workflow):
    """
    This class generate all input files that are needed for
    CHiCAGO to run
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

        logger.info("Generating CHiCAGO input file .rmap")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, input_metadata, output_files):
        """
        This is the main function that run the tools to create .rmap and
        .baitmap.

        Parameters:
        ----------
        input_files: dict
            genome: str
                    Ref genome used in the experiment.

        input_metadata: dict
            input metadata

        output_files: dict
            out_dir_makeRmap: str
                    path to the output diretory
            out_prefix_makeRmap: str
                    prefix for the output file .rmap
            Rtree_files: str
                    Name of the Rtree files

        Returns:
        --------
        bool
        output_metadata_Baitmap : dict
            metadata for both rmap and baitmap
            files
        """

        makeRmap_caller = makeRmapFile(self.configuration)
        output_files_makeRmap, output_metadata_makeRmap = makeRmap_caller.run(
            {
                "genome" : input_files["genome"]
            },
            {
                "genome_digest" : input_metadata["genome_digest"]
            },
            {
                "out_dir_makeRmap" : output_files["out_dir_makeRmap"],
                "out_prefix_makeRmap" : output_files["out_prefix_makeRmap"],
                "Rtree_files" : output_files["Rtree_files"]
            }
        )

        if os.path.getsize(output_files["out_dir_makeRmap"] +
            output_files["out_prefix_makeRmap"] + ".rmap") > 0:
            pass
        else:
            logger.fatal("generate_CHiCAGO_rmap failed to generate .rmap file")
            return False

        return output_files_makeRmap, output_metadata_makeRmap

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
    results = app.launch(generate_CHiCAGO_rmap,
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
        description="Pipeline to generate .rmap file")

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
