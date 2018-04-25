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

from tool.makeDesignFiles import makeDesignFilesTool

#####################################################

class makeDesign(Workflow):
    """
    This class generates the Design files and chinput files,
    imput for CHiCAGO. Starting from rmap and baitmap and capture
    HiC BAM files.
    """

    def __init__(self, configuration=None):
        """
        Initiate the class

        Parameters
        -----------
        Configuration: dict
         dictionary with parameters for different tools from the class
         indicating how to run each of them.
        """

        logger.info("Generating CHiCAGO input Design files")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, input_metadata, output_files):
        """
        Main function to run the tools, MakeDesignFiles_Tool.py and
        bam2chicago_Tool.py

        Parameters
        ----------
        input_files: dict
            designDir: path to the folder with .rmap and .baitmap files
            rmapFile: path to the .rmap file
            baitmapFile: path to the .baitmap file
            bamFile: path to the capture HiC bamfiles

        input_metadata: dict
            input metadata

        output_files: dict
            outPrefixDesign : Path and name of the output prefix,
                recommend to be the same as rmap and baitmap files.
            sample_name: Path and name of the .chinput file

        Returns
        -------
        bool
        output_metadata
        """

        makeDesign_caller = makeDesignFilesTool(self.configuration)
        makeDesgin_out, makeDesign_outMeta = makeDesign_caller.run(
            {
                "designDir" : input_files["designDir"]
            },
            {
                ".rmap" : input_metadata[".rmap"],
                ".baitmap" : input_metadata[".baitmap"]
            },
            {
                "outPrefixDesign" : output_files["outPrefixDesign"]
            }
        )

        if os.path.isfile(output_files["outPrefixDesign"] + ".nbpb") is True:
            pass
        else:
            logger.fatal("processmakeDesign_chinput failed to" +
                         "generate design files")
            return False

        return makeDesgin_out, makeDesign_outMeta

#############################################################

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function

    This function lauch the app using the configuration written
    in two json files:
    """
    #1.Instantiate and launch the app
    print("Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    results = app.launch(makeDesign_chinput,
                         config,
                         in_metadata,
                         out_metadata)

    #2. The App has finished
    print("2. Execution finished: see " +out_metadata)
    print(results)

    return results

#########################################################

if __name__ == "__main__":

    #sert up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Pipeline to generate Design and .chinput files")

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
    IN_METADATA = ARGS.IN_METADATA
    OUT_METADATA = ARGS.OUT_METADATA
    LOCAL = ARGS.local
    if LOCAL:
        import sys
        sys._run_from_cmdl = True

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
