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



from re import compile
from warnings import warn

####################################################

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

        print("makeRmapFile initialising")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)



    def iupac2regex(self, restring):
        """
        Convert target sites with IUPAC nomenclature to regex pattern
        """
        restring = restring.replace('R', '[AG]')
        restring = restring.replace('Y', '[CT]')
        restring = restring.replace('M', '[AC]')
        restring = restring.replace('K', '[GT]')
        restring = restring.replace('S', '[CG]')
        restring = restring.replace('W', '[AT]')
        restring = restring.replace('H', '[ACT]')
        restring = restring.replace('B', '[CGT]')
        restring = restring.replace('V', '[ACG]')
        restring = restring.replace('D', '[AGT]')
        restring = restring.replace('N', '[ATGC]')
        return restring

    def genome_to_dict(self, genome_fa):
        """
        This function takes a genome file in fasta format
        and converts it into a dictionary

        Parameters:
        ----------
        genome_fa : str
            entire path to the genome in fasta file

        Returns:
        --------
            dict
            genome in a dict form, key as chromosomes and
            values str with sequences
        """
        genome_dict = {}
        sequence = ""

        with open(genome_fa, "r") as file:
            for line in file:
                line = line.rstrip()
                if line[0] is ">":
                    if sequence is "":
                        genome_dict[int(line[4:])] = []
                        chromo = int(line[4:])
                        continue
                    else:
                        genome_dict[chromo] = sequence
                        chromo = int(line[4:])
                        sequence = ""
                        continue

                sequence += line
            #Ad last chromosome
            genome_dict[chromo] = sequence

        return(genome_dict)


    def map_re_sites(self, enzyme_name, genome_seq, output_dir, output_prefix, frag_chunk=100000, verbose=False):
        """
        map all restriction enzyme (RE) sites of a given enzyme in a genome.
        Position of a RE site is defined as the genomic coordinate of the first
        nucleotide before the first cut (genomic coordinate starts at 1).
        In the case of HindIII the genomic coordinate is this one:
        123456 789...
             |
             v
        -----A|AGCT T--------------
        -----T TCGA|A--------------
        In this example the coordinate of the RE site would be 6.

        Parameters:
        -----------
        enzyme_name: dict
            name of the enzyme to map (upper/lower case are
            important) as key and value the target sequence
            with a pipe where the enzyme cuts.

        genome_seq: dict
            a dictionary containing the genomic sequence by
            chromosome
        frag_chunk: in order to optimize the search for nearby RE
            sites, each chromosome is splitted into chunks.

        Return
        ------
            bool
            Fragments

        """

        enzymes = enzyme_name

        genome_seq = self.genome_to_dict(genome_seq)

        # we match the full cut-site but report the position after the cut site
        # (third group of the regexp)
        restring = ('%s') % ('|'.join(['(?<=%s(?=%s))' % tuple(enzymes[n].split('|'))
                                     for n in enzymes]))
        # IUPAC conventions
        restring = self.iupac2regex(restring)

        enz_pattern = compile(restring)

        frags = {}
        count = 0
        for crm in genome_seq:
            seq = genome_seq[crm]
            frags[crm] = dict([(i, []) for i in xrange(len(seq) / frag_chunk + 1)])
            frags[crm][0] = [1]
            for match in enz_pattern.finditer(seq):
                pos = match.end()
                frags[crm][pos / frag_chunk].append(pos)
                count += 1
            # at the end of last chunk we add the chromosome length
            frags[crm][len(seq) / frag_chunk].append(len(seq))
            # now we need to assign as first RE site of a fragment the last RE site
            # of previsou fragment, and as last RE site, the first RE site of the
            # next fragment.
            for i in xrange(len(seq) / frag_chunk + 1):
                try:
                    try:
                        frags[crm][i].insert(0, frags[crm][i - 1][-2])
                    except IndexError:
                        # in case there was no RE site in previous fragment
                        frags[crm][i].insert(0, frags[crm][i - 1][-1])
                except KeyError:
                    # it is the very first chunk
                    pass
                plus = 1
                while True:
                    try:
                        frags[crm][i].append(frags[crm][i + plus][0])
                        break
                    except IndexError:
                        # no RE site in this fragment, get "next RE site" from next
                        plus += 1
                    except KeyError:
                        # end of the chromosome
                        break
        if verbose:
            print ('Found' + count + 'RE sites')
        #return frags
        with open(output_dir + output_prefix, "w") as out:
            print(frags, file = out)


    def map_re_sites_nochunk(self, enzyme_name, genome_seq,
        output_dir, output_prefix, verbose=False):
        """
        map all restriction enzyme (RE) sites of a given enzyme in a genome.
        Position of a RE site is defined as the genomic coordinate of the first
        nucleotide before the first cut (genomic coordinate starts at 1).
        In the case of HindIII the genomic coordinate is this one:
        123456 789...
             |
             v
        -----A|AGCT T--------------
        -----T TCGA|A--------------
        In this example the coordinate of the RE site would be 6.


        Parameters:
        -----------
        enzyme_name: dict
            name of the enzyme to map (upper/lower case are
            important) as key and value the target sequence
            with a pipe where the enzyme cuts

        genome_seq: dict
            genome in fasta format

        frag_chunk: in order to optimize the search for nearby RE
            sites, each chromosome is splitted into chunks.

        Return
        ------
            bool
            Fragments
        """
        enzymes = enzyme_name

        #the genome should be in a dictionary
        genome_seq = self.genome_to_dict(genome_seq)

        # we match the full cut-site but report the position after the cut site
        # (third group of the regexp)
        restring = ('%s') % ('|'.join(['(?<=%s(?=%s))' % tuple(enzymes[n].split('|'))
                                    for n in enzymes]))
        # IUPAC conventions
        restring = self.iupac2regex(restring)

        enz_pattern = compile(restring)

        frags = {}
        count = 0
        for crm in genome_seq:
            seq = genome_seq[crm]
            frags[crm] = [1]
            for match in enz_pattern.finditer(seq):
                pos = match.end()
                frags[crm].append(pos)
                count += 1
            # at the end of last chunk we add the chromosome length
            frags[crm].append(len(seq))
        if verbose:
            print ('Found %d RE sites' % count)

        with open(output_dir + output_prefix, "w") as out:
            print(frags, file = out)

    #def from_frag_to_rmap(self, frags, output_dir, output_prefix):



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

        results = self.map_re_sites_nochunk(
            input_files["RE"],
            input_files["genome"],
            output_files["output_dir"],
            output_files["output_prefix"]
            )

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
    "genome" : "../genomes/GRCh38/GRCh38.fa",
    "RE" : { "HindIII" : 'A|AGCTT'}
}


metadata = {"genome" : Metadata(
    "hg38", "fasta", "../genome_mm10/mm10.fa", None, "HindIII", 9606),
}


output_files = {
    "output_dir" : "../tests/data/test_makeRmap/",
    "output_prefix" : "restriction_enzyme_test.txt"
}

test = makeRmapFile()

test.run(input_files, metadata, output_files)

#print(test.run(input_files, metadata, output_files))
