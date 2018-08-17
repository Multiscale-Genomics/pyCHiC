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

import argparse

from basic_modules.workflow import Workflow
from utils import logger

from CHiC.tool.rmap_tool import makeRmapFile

################################################

class process_rmap(Workflow):
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

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that run the tools to create .rmap and
        .baitmap.

        Parameters:
        ----------
        input_files: dict
            genome_fa: str
                    Ref genome used in the experiment.

        metadata: dict
            input metadata

        output_files: dict
            out_dir_rmap: str
                    path to the output diretory
            out_prefix_rmap: str
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
                    "Rtree_file_idx" : output_files["Rtree_file_idx"],
                    "chr_handler" : output_files["chr_handler"]
                }
            )



            logger.info(".rmap file generated succesfully")

        except IOError:
            logger.fatal("rmap_tool failed to generate .rmap file")


        return output_files_rmap, output_metadata_rmap


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
    results = app.launch(process_rmap,
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
        description="Pipeline to generate .rmap file")

    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument(
        "--in_metadata", help="Location of metadata file")
    PARSER.add_argument(
        "--out_metadata", help="Location of output metadata file")
    PARSER.add_argument(
        "--local", action="store_const", const=True, default=False)

    #Get matching Parametersmeters from the command line
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
