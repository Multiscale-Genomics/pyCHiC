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
import subprocess
import sys

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports


from basic_modules.tool import Tool
from basic_modules.metadata import Metadata
###############################################################

class makeRmapFile(Tool):
    """
    Tool for digest the genome with one RE. Wrapper of hicup_digester
    """

    def __init__(self, configuration=None):
        """
        initialising the function

        Parameters:
        -----------
        configuration: dict
         dictionary containing all the arguments and parameters
         to run the tool
        """

        print("bam2chicago initialising")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def hicup_digester(self, genome, arguments, output_dir, output_prefix):
        """
        This function runs the Perl script hicup_digester

        Parameters
        ----------
        genome : str
            path to the genome to be digested
        args: list
            Select the args that the user define

        Output:
        --------
        output: str
            path to the output file
        bool
        """


        args = ["hicup_digester", genome,
            "--outdir", output_dir, "--genome", output_prefix]

        args = args + arguments

        print(args)
        logger.info("hicup_digester CMD: " + " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr= subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            if "digester_re2" in self.configuration:
                os.path.isfile("Digest_"+ output_dir + output_prefix +
                "_" + self.configuration["digester_re1"].split(",")[1]
                + "_" + self.configuration["digester_re2"].split(",")[1]
                + ".txt")
            else:
                 os.path.isfile("Digest_"+ output_dir + output_prefix +
                "_" + self.configuration["digester_re1"].split(",")[1]
                + "_None.txt")

        except IOError:
            logger.fatal("chicago failed to generate peak file")
            logger.fatal("chicago stdout" + proc_out)
            logger.fatal("chicago stderr" + proc_err)
            return False

        return True

    #def adjust_format_to_rmpa(output_digest, output_formated):
        """
        This function takes the output file from hicup_digester
        and put it in the correct format for rmap chicago input.
        <chr> <start> <end> <numeric ID>
        """


    @staticmethod
    def get_digester_params(params):
        """
        Function that select the parameters that have been
        selected by the user and pass them to run hicup_digester

        Parameters
        ----------
        params: dict
            --re1       Restriction enzyme used to digest the genome (the enzyme that
                        forms the ligation junction) e.g. A^GATCT,BglII.  Some Hi-C protocols may use two enzymes at this stage.  To specify two enzymes: -1 A^GATCT,BglII:A^AGCTT,HindIII.
            --re2       To specify a restriction enzyme instead of sonication to
                        shorten di-tags. This restriction site does NOT form a Hi-C ligation junction. 2 .g. AG^CT,AluI. Typically the sonication protocol is followed.
            --config    Specify the name of the optional configuration file
            --genome    Name of the genome to be digested (not the path to the genome
                        file or files, but the genome name to include in the output file)
            --help      Print program help and exit
            --outdir    Specify the directory to which the output files should be
                        written
            --quiet     Suppress all progress reports
            --version   Print the program version and exit
            --zip       Print the results to a gzip file
        """

        command_params = []

        command_parameters = {
            "digester_re1" : ["--re1", True],
            "digester_re2" : ["--re2", True],
            "digester_quite" :["--quiet", False],
            "digester_version" : ["--version", False],
            "digester_zip" : ["--zip", False]
            }

        for arg in params:
            if arg in command_parameters:
                if command_parameters[arg][1] is True:
                    command_params += [command_parameters[arg][0], params[arg]]
                else:
                    command_params += [command_parameters[arg][0]]
        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        This function run the tool

        Parameters
        ----------

        input_files: dict
            genome in fasta file
        input_metadata: dict
            input metadata

        Returns
        -------

        output_files: dict
            name and path to the output file
        output_metadata: dict
            lest of matching metadata
        """

        command_params = self.get_digester_params(self.configuration)

        logger.info("hicup_digester command parameters "+
         " ".join(command_params))

        results = self.hicup_digester(input_files["genome"],
            command_params, output_files["output_dir"], output_files["output_prefix"])

        results = compss_wait_on(results)

        output_metadata = {
        }

        return(results, output_metadata)

#from the command line  ../hicup_digester --re1 A^AGCTT,HindIII --genome caquita_digestiva --outdir . /Users/pacera/developing/C-HiC/genome_mm10/mm10.fa


input_files = {
    "genome" : "../genome_mm10/mm10.fa"
}


metadata = {

}

config = {
     "digester_re1" : "A^AGCTT,HindIII"
}

output_files = {
    "output_dir" : ".",
    "output_prefix" : "test_digest"
}

test = makeRmapFile(config)

test.run(input_files, metadata, output_files)




