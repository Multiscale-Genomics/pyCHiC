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
from shutil import copy
from rtree import index

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
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata
from tool.bowtie_aligner import bowtie2AlignerTool

##################################################


class makeBaitmapTool(Tool):
    """
    This tool use probe capture sequences as input.
    Then using bwa can tell which baits correspond
    to the probes
    """
    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration

        Parameters
        ----------
        configuration: dict
            parameters to run the tool
        """

        logger.info("initialising makeBaitmapTool")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=list, sam_file=FILE_OUT, out_bam=FILE_IN, rtree_dat=FILE_IN, rtree_idx=FILE_IN,
          rtree_prefix=IN, chr_handler=FILE_IN)
    def sam_to_baitmap(self, sam_file, out_bam, rtree_dat, rtree_idx, rtree_prefix,
                       chr_handler):  # pylint: disable=no-self-use
        """
        This function take the sam file, output of bwa
        and the Rtree_files, and output a baitmap file

        Parameters
        ----------
        sam_file : str
            path to output file from bwa_for_probes
            complete path to .rmap file
        """

        args = ["samtools", "view", "-h", "-o", sam_file, out_bam]

        logger.info("samtools args: " + ' '.join(args))

        try:
            with open(sam_file, "w") as f_out:
                process = subprocess.Popen(
                    ' '.join(args),
                    shell=True,
                    stdout=f_out, stderr=f_out
                    )
            process.wait()

        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, args))
            return False

        chr_dict = {}
        with open(chr_handler, "r") as chr_file:
            for line in chr_file:
                line_hdl = line.rstrip().split("\t")
                chr_dict[line_hdl[1]] = int(line_hdl[0])

        copy(rtree_idx, rtree_prefix+".idx")
        copy(rtree_dat, rtree_prefix+".dat")

        idx = index.Rtree(rtree_prefix)

        baitmap = []

        features = []

        with open(sam_file, "r") as file_in:
            for line in file_in:
                if line[0] != "@":
                    line = line.rstrip().split("\t")

                    try:
                        crm = chr_dict[line[2]]
                    except IOError:
                        continue

                    features.append(line[0])

                    srt_pos = int(line[3])
                    end_pos = srt_pos + int(len(line[9]))

                    id_object = idx.intersection(
                        (srt_pos, crm, end_pos, crm),
                        objects=True)

                    hits = [[i.id, i.bbox] for i in id_object]

                    if len(hits) > 1:
                        logger.warning("probe map to two RE fragmnets, " +
                                       " ".join(line) + " start pos" + str(srt_pos) +
                                       " end pos" + str(end_pos))

                    elif not hits:
                        logger.warn("Sequence does not" +
                                    "match with any RE fragment, " +
                                    " ".join(line))
                        continue

                    else:
                        fragment_coord = [
                            int(hits[0][1][1]),
                            int(hits[0][1][0]),
                            int(hits[0][1][2]),
                            hits[0][0]
                            ]

                        if fragment_coord not in baitmap:
                            baitmap.append(fragment_coord)

        os.remove(rtree_prefix+".dat")
        os.remove(rtree_prefix+".idx")

        return baitmap, features

    @task(returns=bool, features=IN, baitmap_list=IN,
          out_baitmap=FILE_OUT, chr_handler=FILE_IN)
    def create_baitmap(self, features, # pylint: disable=no-self-use
                       baitmap_list,
                       out_baitmap,
                       chr_handler):
        """
        This function takes a list with RE fragments that
        correspond to baits and print it to a file

        Parameters
        ----------
        baitmap_list: list
            lsit with all the RE fragments corresponding
            to baits
        out: str
            entire pat and name of the .baitmap file
        """
        # print(out_baitmap)

        chr_dict = {}
        with open(chr_handler, "r") as chr_file:
            for line in chr_file:
                line_hdl = line.rstrip().split("\t")
                chr_dict[int(line_hdl[0])] = line_hdl[1]

        with open(out_baitmap, "a") as file_out:
            for i in enumerate(baitmap_list):
                file_out.write("{}\t{}\t{}\t{}\t{}\n".format(
                    chr_dict[i[1][0]],
                    i[1][1],
                    i[1][2],
                    i[1][3],
                    features[i[0]]
                    ))

        if os.path.getsize(out_baitmap) > 0:
            return True

        logger.fatal("baitmap file not generated")
        return False

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to produce a .baitmap file, starting from rtree files,
        indexed genome and probes.

        Parameters
        ----------
        input_file : dict
            a dict of absolute path names of the input data elements,
            associated with their role;
        input_metadata : dict
            a dict of metadatas for each of the input data elements,
            associated with their role;
        output_files : dict
            a dict of absolute path names of the output data elements,
            associated with their role.

        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        """

        if "genome_fa_public" in input_files:
            input_files["genome_fa"] = input_files.pop("genome_fa_public")
            input_metadata["genome_fa"] = input_metadata.pop("genome_fa_public")

            input_files["bowtie_gen_idx"] = input_files.pop("bowtie_gen_idx_public")
            input_metadata["bowtie_gen_idx"] = input_metadata.pop("bowtie_gen_idx_public")


        if os.path.isfile("bam.bai"):
            os.remove("bam.bai")

        out_bam = "tests/data/test_baitmap/baits.bam"
        rtree_file_dat = "tests/data/test_rmap/rtree_file.dat"
        rtree_file_idx = "tests/data/test_rmap/rtree_file.idx"
        chr_handler = "tests/data/test_baitmap/chr_handler.txt"
        RMAP = "tests/data/test_run_chicago/test.rmap"
        out_baitmap = "tests/data/test_run_chicago/test.baitmap"
        bait_sam = "tests/data/test_baitmap/baits.sam"

        re_meta = {
            self.configuration["chic_RE_name"]: self.configuration["chic_RE_sequence"],
            "bowtie2_fasta_input" : "True"}

        input_bwa = {
            "genome": input_files["genome_fa"],
            "index": input_files["bowtie_gen_idx"],
            "loc": input_files["probes_fa"]
        }

        output_bwa = {
            "output": out_bam,
            "bai" : "bam.bai"
        }
        metadata_bwa = {
            "genome": input_metadata["genome_fa"],
            "index": input_metadata["bowtie_gen_idx"],
            "loc": input_metadata["probes_fa"]
        }

        bowtie2_t = bowtie2AlignerTool(self.configuration)
        bowtie_files, bowtie_meta = bowtie2_t.run(input_bwa, metadata_bwa, output_bwa)

        # bwa_meta = compss_wait_on(bwa_meta)


        prefix_rtree = "rtree_file"

        baitmap_list, features = self.sam_to_baitmap(
            bait_sam,
            bowtie_files["bam"],
            rtree_file_dat,
            rtree_file_idx,
            prefix_rtree,
            chr_handler)

        results = self.create_baitmap(
            features,
            baitmap_list,
            out_baitmap,
            chr_handler)

        #results = compss_wait_on(results)

        output_metadata = {
            }


        return output_files, output_metadata
