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

from CHiC.tool.makeBaitmap import makeBaitmapTool

################################################


class process_baitmap(Workflow):
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

        logger.info("Generating CHiCAGO input file .baitmap")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that run the tools to create
        .baitmap.

        Parameters
        ----------
        input_files: dict
            genome: str
                Indexed Ref genome used used in the experiment.

        metadata: dict
            input metadata

        output_files: dict
            Rtree_files: str
                    Name of the Rtree files
            bait_sam: str
                whole path of SAM file, generated to locate baits
            out_baitmap: str
                whole path for the .baitmap file

        Returns
        -------
        output_files_Baitmap : bool
        output_metadata_Baitmap : dict
            metadata for baitmap
            files
        """
        try:
            baitmap_caller = makeBaitmapTool(self.configuration)
            output_files_baitmap, output_metadata_baitmap = baitmap_caller.run(
                {
                    "genome_idx": input_files["genome_idx"],
                    "probes_fa": input_files["probes_fa"],
                    "genome_fa": input_files["genome_fa"],
                },
                {
                    "genome_fa": metadata["genome_fa"],
                    "probes_fa": metadata["probes_fa"],
                    "genome_idx": metadata["genome_idx"]
                },
                {

                }
            )

            logger.info(".baitmap file generated succesfully")

        except IOError:
            logger.fatal("generate_CHiCAGO_baitmap failed to generate .baitmap file")

        return output_files_baitmap, output_metadata_baitmap

#############################################################


def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function

    This function launch the app using the configuration written
    in two json files: config_process_rmapBaitmap.json  and
    input_process_rmapBaitmap.json
    """
    # Instantiate and lauch the app
    print("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    results = app.launch(process_baitmap,
                         config,
                         in_metadata,
                         out_metadata)
    # 2. The App has finished
    print("2. Execution finished: see " + out_metadata)
    print(results)

    return results

#########################################################################


if __name__ == "__main__":

    # set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Pipeline to generate .baitmap file")

    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument(
        "--in_metadata", help="Location of metadata file")
    PARSER.add_argument(
        "--out_metadata", help="Location of output metadata file")
    PARSER.add_argument(
        "--local", action="store_const", const=True, default=False)

    # Get matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata
    LOCAL = ARGS.local

    if LOCAL:
        import sys
        sys._run_from_cmdl = True  # pylint: disable=protected-access

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
