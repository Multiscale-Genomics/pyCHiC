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
from re import compile

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

    @staticmethod

    def map_re_sites_nochunk(self, enzyme_name, genome_seq, verbose=False):
    """
    map all restriction enzyme (RE) sites of a given enzyme in a genome.
    Position of a RE site is defined as the genomic coordinate of the first
    nucleotide after the first cut (genomic coordinate starts at 1).
    In the case of HindIII the genomic coordinate is this one:
    123456 789...
           |
           v
    -----A|AGCT T--------------
    -----T TCGA|A--------------
    In this example the coordinate of the RE site would be 7.
    :param enzyme_name: name of the enzyme to map (upper/lower case are
       important)
    :param genome_seq: a dictionary containing the genomic sequence by
       chromosome
    """

    if isinstance(enzyme_name, str):
        enzyme_names = [enzyme_name]
    elif isinstance(enzyme_name, list):
        enzyme_names = enzyme_name

    enzymes = {}

    for name in enzyme_names:
        enzymes[name] = RESTRICTION_ENZYMES[name]

    # we match the full cut-site but report the position after the cut site
    # (third group of the regexp)
    restring = ('%s') % ('|'.join(['(?<=%s(?=%s))' % tuple(enzymes[n].split('|'))
                                   for n in enzymes]))
    # IUPAC conventions
    restring = iupac2regex(restring)

    enz_pattern = compile(restring)

    frags = {}
    count = 0
    for crm in genome_seq:
        seq = genome_seq[crm]
        frags[crm] = [1]
        for match in enz_pattern.finditer(seq):
            pos = match.end() + 1
            frags[crm].append(pos)
            count += 1
        # at the end of last chunk we add the chromosome length
        frags[crm].append(len(seq))
    if verbose:
        print 'Found %d RE sites' % count
    return frags





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

        process = subprocess.Popen(args,
            stdout=subprocess.PIPE,
            stderr= subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        out = "".join([f for f in os.listdir(output_dir)
            if f.startswith("Digest_")])

        if os.path.getsize(output_dir + out) > 0:
            pass
        else:
            logger.fatal("hicup_digester failed to generate peak file")
            logger.fatal("hicup_digester stdout" + proc_out)
            logger.fatal("hicup_digester stderr" + proc_err)
            return False

        if self.adjust_format_to_rmap(output_dir, out) is True:
            return True


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
            "rmap": Metadata(
                data_type=input_metadata['genome'].data_type,
                file_type="rmap",
                file_path=output_files["output_dir"],
                sources=[
                    input_metadata["genome"].file_path,
                ],
                taxon_id=input_metadata["genome"].taxon_id,
                meta_data={
                    "RE" : input_metadata["genome"].meta_data,
                    "tool": "makeRmapFileTool"
                }
            )

        }

        return(results, output_metadata)


input_files = {
    "genome" : "../genomes/GRCh38/GRCh38.fa"
}


metadata = {"genome" : Metadata(
    "hg38", "fasta", "../genome_mm10/mm10.fa", None, "HindIII", 9606),
}

config = {
     "digester_re1" : "A^AGCTT,HindIII"
}

output_files = {
    "output_dir" : "../tests/data/test_makeRmap/",
    "output_prefix" : "hicup_digester_GRCh38"
}

test = makeRmapFile(config)

test.run(input_files, metadata, output_files)




