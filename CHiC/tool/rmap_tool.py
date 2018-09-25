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
import sys
import re
from shutil import move

from rtree import index

from utils import logger
from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

re.compile("pattern")

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

####################################################

class makeRmapFile(Tool):

    """
    Tool for digest the genome with one renzime. Wrapper of hicup_digester
    """


    def __init__(self, configuration=None):
        """
        initialising the function

        Parameters
        ----------
        configuration: dict
            dictionary containing all the arguments and parameters
            to run the tool
        """

        print("makeRmapFile initialising")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def iupac2regex(restring):
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

    @staticmethod
    def genome_to_dict(genome_fa):
        """
        This function takes a genome file in fasta format
        and converts it into a dictionary

        Parameters
        ----------
        genome_fa : str
            entire path to the genome in fasta file

        Returns
        -------
            dict
            genome in a dict form, key as chromosomes and
            values str with sequences
        """
        logger.info("converting fasta genome into a dictionary")
        genome_dict = {}
        sequence = ""

        chromo_fake = 1
        chromo_dict = {}

        with open(genome_fa, "r") as file_handle:
            for line in file_handle:
                line = line.rstrip()
                if line[0] == ">":
                    if not sequence:
                        chromo_dict[chromo_fake] = line[1:]
                        continue
                    else:
                        genome_dict[chromo_fake] = sequence
                        chromo_fake += 1
                        chromo_dict[chromo_fake] = line[1:]
                        sequence = ""
                        continue

                sequence += line.upper()
            #Ad last chromosome
            genome_dict[chromo_fake] = sequence

        return genome_dict, chromo_dict

    def map_re_sites2(self, enzyme_name, genome_fa):
        """
        map all restriction enzyme (renzime) sites of a given enzyme in a genome.
        Position of a renzime site is defined as the genomic coordinate of the first
        nucleotide before the first cut (genomic coordinate starts at 1).
        In the case of HindIII the genomic coordinate is this one:

.. code-block:: none

        123456 789...
             |
             v
        -----A|AGCT T--------------
        -----T TCGA|A--------------
        In this example the coordinate of the renzime site would be 6.

        Parameters
        ----------
        enzyme_name: dict
            name of the enzyme to map (upper/lower case are
            important) as key and value the target sequence
            with a pipe where the enzyme cuts

        genome_fa: str
            genome in fasta format

        Returns
        -------
            list
        """
        try:
            enzymes = enzyme_name

            #the genome should be in a dictionary
            genome_seq, chromo_dict = self.genome_to_dict(genome_fa)

            # we match the full cut-site but report the position after the cut site
            # (third group of the regexp)
            restring = ('%s') % (
                '|'.join(
                    [
                        '(?<=%s(?=%s))' % tuple(enzymes[n].split('|'))
                        for n in enzymes]
                    )
                )

            # IUPAC conventions
            restring = self.iupac2regex(restring)

            enz_pattern = re.compile(restring)

            frags = {}
            count = 0

            logger.info("searching renzime sites")

            for crm in genome_seq:
                seq = genome_seq[crm]
                frags[crm] = []
                for match in enz_pattern.finditer(seq):
                    pos = match.end()
                    frags[crm].append(pos)
                    count += 1
                # at the end of last chunk we add the chromosome length
                frags[crm].append(len(seq))

            return frags, chromo_dict

        except IOError:
            logger.fatal("map_re_sites2 function from rmap_tool failed =(")

    @task(returns=bool, enzyme_name=IN, genome_fa=FILE_IN,
          rtree=IN, rtree_dat=FILE_OUT, rtree_idx=FILE_OUT,
          RMAP=FILE_OUT, chr_handler=FILE_OUT)
    def from_frag_to_rmap(self, enzyme_name, genome_fa,
                          rtree, rtree_dat, rtree_idx, RMAP, chr_handler):
        """
        This function takes the fragment output from digestion and
        convert them into rmap files.

        It also save the renzime sites positions and ID into a file using Rtree
        python module. This file will be used by makeBatmap.py to generate
        .batmap file using spatial indexing

        Parameters
        ----------
        enzyme_name: str
            described in map_re_sites2
        genome_fa: str
            full path to genome FASTA format
        frags : dict
            dict containing chromosomes as keys and
            renzime sites as values
        out_dir_rmap: str
            path to the output directory
        out_prefix_rmap: str
            name of the output file.
        """
        #include creation folders
        frags, chromo_dict = self.map_re_sites2(enzyme_name, genome_fa)

        logger.info("coverting renzime fragments into rmap file")

        try:
            idx = index.Rtree(rtree)
        except AttributeError:
            logger.info("index failed =(")

        with open(RMAP, "w") as out:
            counter_id = 0
            for crm in frags:
                counter = 0
                for re_site in frags[crm]:
                    counter_id += 1
                    counter += 1
                    if counter == 1:
                        out.write("{}\t{}\t{}\t{}\n".format(chromo_dict[crm],
                                                            1,
                                                            re_site,
                                                            counter_id),
                                 )
                        idx.insert(counter_id, (1, crm, re_site, crm))
                    else:
                        out.write("{}\t{}\t{}\t{}\n".format(chromo_dict[crm],
                                                            prev_re_site+1, # pylint: disable=used-before-assignment
                                                            re_site,
                                                            counter_id),
                                 )
                        idx.insert(counter_id, (prev_re_site+1, crm, re_site, crm))

                    prev_re_site = re_site

        idx.close()

        with open(chr_handler, "w") as chr_file:
            for chr_fake, chr_real in chromo_dict.items():
                chr_file.write("{}\t{}\n".format(chr_fake, chr_real))

        try:
            move(rtree+".dat", rtree_dat)
            move(rtree+".idx", rtree_idx)
            return True

        except IOError:
            logger.fatal("makeRmap_Tool.py failed to generate .rmap file")
            return False


    def run(self, input_files, metadata, output_files):
        """
        This function run the tool

        Parameters
        ----------

        input_files: dict
            genome in fasta file
        metadata: dict
            input metadata

        Returns
        -------

        output_files: dict
            name and path to the output file
        output_metadata: dict
            lest of matching metadata
        """
        out_dir_rmap = "/".join(output_files["RMAP"].split("/")[:-1])
        out_dir_rtree = "/".join(output_files["Rtree_file_dat"].split("/")[:-1])

        if not os.path.isdir(out_dir_rmap):
            os.mkdir(out_dir_rmap)

        if not os.path.isdir(out_dir_rtree):
            os.mkdir(out_dir_rtree)


        if "".join(output_files["Rtree_file_dat"].split(".")[:-1]) != \
           "".join(output_files["Rtree_file_idx"].split(".")[:-1]):
            logger.fatal("Rtree_file_dat and Rtree_file_idx"
                         "should have the same prefix name")

        rtree = "".join(output_files["Rtree_file_dat"].split(".")[:-1])

        rtree = "".join(rtree.split("/")[-1])

        results = self.from_frag_to_rmap(
            self.configuration["renzime"],
            input_files["genome_fa"],
            rtree,
            output_files["Rtree_file_dat"],
            output_files["Rtree_file_idx"],
            output_files["RMAP"],
            output_files["chr_handler"]
        )

        #results = compss_wait_on(results)

        output_metadata = {
            "RMAP": Metadata(
                data_type=metadata['genome_fa'].data_type,
                file_type="rmap",
                file_path=output_files["RMAP"],
                sources=[
                    metadata["genome_fa"].file_path
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "renzime" : self.configuration["renzime"],
                    "tool": "rmap_tool"
                }
            ),
            "Rtree_file_dat": Metadata(
                data_type="Rtree_file_dat",
                file_type="Rtree_file_dat",
                file_path=output_files["Rtree_file_dat"],
                sources=[
                    metadata["genome_fa"].file_path
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "renzime" : self.configuration["renzime"],
                    "tool": "rmap_tool"
                }
            ),
            "Rtree_file_idx": Metadata(
                data_type="Rtree_file_idx",
                file_type="Rtree_file_idx",
                file_path=output_files["Rtree_file_dat"],
                sources=[
                    metadata["genome_fa"].file_path
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "renzime" : self.configuration["renzime"],
                    "tool": "rmap_tool"
                }
            ),
            "chr_handler": Metadata(
                data_type="data_CHi-C",
                file_type="TXT",
                file_path=output_files["chr_handler"],
                sources=[
                    metadata["genome_fa"].file_path
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "renzime" : self.configuration["renzime"]
                }
            )


        }

        return output_files, output_metadata
