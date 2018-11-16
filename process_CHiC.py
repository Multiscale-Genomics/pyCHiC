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
import os
from shutil import move
from shutil import rmtree

from basic_modules.workflow import Workflow
from utils import logger


from CHiC.tool.rmap_tool import makeRmapFile
from CHiC.tool.makeBaitmap import makeBaitmapTool
from CHiC.tool.makeDesignFiles import makeDesignFilesTool
from CHiC.tool.hicup_tool import hicup
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

        Parameters
        ----------
        input_files: dict
            fastq1: str
            fastq2: str
            genome_fa: str
            genome in fasta format

        input_metadata: dict
            input metadata

        output_files: dict

        Returns
        -------
        output_files
        output_metadata

        """
        self.configuration["bowtie2_fasta_input"] = "True"

        if "genome_fa_public" in input_files:
            input_files["genome_fa"] = input_files.pop("genome_fa_public")
            metadata["genome_fa"] = metadata.pop("genome_fa_public")

            input_files["bowtie_gen_idx"] = input_files.pop("bowtie_gen_idx_public")
            metadata["bowtie_gen_idx"] = metadata.pop("bowtie_gen_idx_public")


        #call hicup
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

            logger.info("hicup runned succesfully =)")

        except IOError:
            logger.fatal("hicup failed to run succesfully =(")

        #call ramp
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
                }
            )
            logger.info(".rmap file generated succesfully")

        except IOError:
            logger.fatal("rmap_tool failed to generate .rmap file")

        #produce Baitmap file
        try:
            baitmap_caller = makeBaitmapTool(self.configuration)
            output_files_baitmap, output_metadata_baitmap = baitmap_caller.run(
                {
                    "bowtie_gen_idx" : input_files["bowtie_gen_idx"],
                    "probes_fa": input_files["probes_fa"],
                    "genome_fa" : input_files["genome_fa"],
                },
                {
                    "genome_fa" : metadata["genome_fa"],
                    "probes_fa" : metadata["probes_fa"],
                    "bowtie_gen_idx": metadata["bowtie_gen_idx"]
                },
                {
                }
            )

            logger.info(".baitmap file generated succesfully")

        except IOError:
            logger.fatal("generate_CHiCAGO_baitmap failed to generate .baitmap file")


        try:
            design_caller = makeDesignFilesTool(self.configuration)
            design_out, design_meta = design_caller.run(
                {

                },
                {
                },
                {
                }
            )

            logger.info("design files succesfully generated =)")

        except IOError:
            logger.fatal("process_makeDesign failed to" +
                         "generate design files")

        try:
            bam2chicago_caller = bam2chicagoTool(self.configuration)
            output_files_bam2chicago, output_metadata_bam2chicago = bam2chicago_caller.run(
               input_files, metadata, output_files
            )

            logger.info("bam2chicago_tool succesfully generate chinput files =)")

        except IOError:
            logger.fatal("process_bam2chicago failed to generate .chinput files")

        try:
            chicago_caller = ChicagoTool(self.configuration)

            output_files_chicago, output_metadata_chicago = chicago_caller.run(
                input_files, metadata, output_files)

        except IOError:
            logger.info("chicago failed to generate output files =(")

        output_files = {}
        output_metadata = {}

        output_files.update(output_files_rmap)
        output_files.update(output_files_baitmap)
        output_files.update(design_out)
        output_files.update(output_files_hicup)
        output_files.update(output_files_bam2chicago)
        output_files.update(output_files_chicago)

        output_metadata.update(output_metadata_rmap)
        output_metadata.update(output_metadata_baitmap)
        output_metadata.update(design_meta)
        output_metadata.update(output_metadata_hicup)
        output_metadata.update(output_metadata_bam2chicago)
        output_metadata.update(output_metadata_chicago)


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
