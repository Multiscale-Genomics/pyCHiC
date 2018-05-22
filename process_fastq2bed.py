#!/usr/bin/env python

"""
    See the NOTICE file distributed with this work for additional information
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
import os

from basic_modules.workflow import Workflow
from utils import logger

from tool.fastq2bed import Fastq2bed

#################################################

class process_Fastq2bed(Workflow):
    """
    This class convert pair fastq reads to bed files
    that will feed process_bed2chicagobamWrap.py
    """

    def __init__(self, configuration=None):
        """
        Initialising the class

        Parameter
        ---------
        configuration: dict
            dictionary with all the parameters tu run
            the functions
        """

        logger.info("initialising process_fastq2bed")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, input_metadata, output_files):
        """
        this function runs all functions from fastq2bed.py tool

        Parameters
        ----------
        input_files: dict
            fastq1: str
                path to fastq reads 1
            fastq2: str
                path to fastq reads 2
            RE: str
                restriction enzyme, used to digest the genome.
                format is the name of the eznyme case sensitive
                example HindIII
            chr: str
                chromosome number to run just one chromosome
            genindex: str
                path to the genome indexed using GEM
            genome_fasta: str
                path to the genome in fasta format

        input_metadata: dic
            metadata from input

        output_files: dict
            wd:str
                path to the working directory
                where there are the output files

        Returns
        -------
        Bool
        output_metadata: dict
            metadata for fastq2bed output
        """
        fastq2bed_caller = Fastq2bed(self.configuration)
        out_fastq2bed, out_meta_fastq2bed = fastq2bed_caller.run(
            {
                "fastq1" : input_files["fastq1"],
                "fastq2" : input_files["fastq2"],
                "RE" : input_files["RE"],
                "chromosome": input_files["chromosome"],
                "genindex": input_files["genindex"],
                "genome_fasta" : input_files["genome_fasta"]
            },{
                "fastq1": input_metadata["fastq1"],
                "fastq2": input_metadata["fastq2"],
                "genome_fasta": input_metadata["genome_fasta"]
            },{
                "wd": output_files["wd"]
            }
        )

        try:
            file_list = os.listdir(output_files["wd"]+"/03_filtered_reads/")
            for file_hdl in file_list:
                if "valid" in file_hdl:
                    assert os.path.getsize(output_files["wd"]+"/03_filtered_reads/"+file_hdl) > 0
        except IOError:
            print("process_fastq2bed didn't generate any output =(")

        return out_fastq2bed, out_meta_fastq2bed

###################################################################

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function

    This function launch the app using the configuration written
    in two json files: config_process_rmapBaitmap.json  and
    input_process_rmapBaitmap.json
    """
    #Instantiate and lauch the app
    print("1.Instantiate and launch the App")
    from app.jsonapp import JASONApp
    app = JASONApp()
    results = app.launch(process_Fastq2bed,
                        config, in_metadata,
                        out_metadata)
    #2. Th2 App has finished
    print("2. Execution is finished: " + out_metadata)
    print(results)

    return results

##################################################################


if __name__ == "__name__":

    #set up the command line parameters
    PARSER = argparse.ArgumentParser(
        description="Pipeline to generate bed (.tsv) files")

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