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

from CHiC.tool.bed2bam import bed2bam

################################################

class process_bed2bam(Workflow):
    """
    This class generate bam file compatible with
    bam2chicago.py from .tsv file ourput of process_fastq2bed.py
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

        logger.info("Initialising process_bed2bam")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that run the tools to create bamfile

        Parameters:
        ----------
        input_files: dict
            bed : str
                path to the bed file


        metadata: dict
            input metadata

        output_files: dict
            bam_out : str
                complete path and prefix of output

        Returns:
        --------
        bool
        output_metadata_Baitmap : dict
            metadata for both rmap and baitmap
            files
        """

        bed2bam_caller = bed2bam(self.configuration)
        output_files_bed2bam, output_metadata_bed2bam = bed2bam_caller.run(
            {
                "bed" : input_files["bed"]
            },
            {
                "bed" : metadata["bed"]
            },
            {
                "bam_out" : output_files["bam_out"],
                "bam_out_sorted" : output_files["bam_out_sorted"]
            }
        )

        if os.path.getsize(output_files["bam_out_sorted"]) > 0:
            pass
        else:
            logger.fatal("process_bed2bam failed to generate BAM file")
            return False

        return output_files_bed2bam, output_metadata_bed2bam

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
    results = app.launch(process_bed2bam,
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
        description="Pipeline to generate BAM files")

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
