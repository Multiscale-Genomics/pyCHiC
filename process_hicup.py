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

from CHiC.tool.hicup_tool import hicup

################################################

class process_hicup(Workflow):
    """
    This class run hicup tool which run hicup, doing the alignment and
    filtering of the reads and convert them into a BAM file.
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

        logger.info("Truncate C-HiC fastq reads")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
        """
        This is the main function that runs

        Parameters
        ----------
        input_files: dict
            fastq1
            fastq2

        metadata: dict

        output_files: dict
            out_dir: str
                directory to write the output

        Returns
        -------
        results: bool
        output_metadata: dict
        """
        if "genome_fa_public" in input_files:
            input_files["genome_fa"] = input_files.pop("genome_fa_public")
            metadata["genome_fa"] = metadata.pop("genome_fa_public")

            input_files["bowtie_gen_idx"] = input_files.pop("bowtie_gen_idx_public")
            metadata["bowtie_gen_idx"] = metadata.pop("bowtie_gen_idx_public")

        try:
            hicup_caller = hicup(self.configuration)
            output_files_hicup, output_metadata_hicup = hicup_caller.run(
                {
                    "genome_fa": input_files["genome_fa"],
                    "fastq1": input_files["fastq1"],
                    "fastq2" : input_files["fastq2"],
                    "bowtie_gen_idx": input_files["bowtie_gen_idx"]
                },
                {
                    "genome_fa": metadata["genome_fa"],
                    "fastq1": metadata["fastq1"],
                    "fastq2": metadata["fastq2"]
                },
                {
                    "hicup_outdir_tar" : output_files["hicup_outdir_tar"]
                }
            )


            #if os.path.isfile(output_files["fastq1_trunc"]) is True:
                #if os.path.getsize(output_files["fastq1_trunc"]) > 0:
            return output_files_hicup, output_metadata_hicup

        except IOError:
            return False

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
    results = app.launch(process_hicup,
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
        description="Pipeline to truncate fastq reads")

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
