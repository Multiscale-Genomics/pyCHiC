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
import shutil
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
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata
from tool.bwa_mem_aligner import bwaAlignerMEMTool

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

        Parameters:
        -----------
        configuration: dict
            parameters to run the tool
        """

        logger.info("initialising makeBaitmapTool")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)


    @task(returns=list, sam_file=FILE_OUT, out_bam=FILE_IN, rtree_dat=FILE_IN, rtree_idx=FILE_IN,
          rtree_prefix=IN)
    def sam_to_baitmap(self, sam_file, out_bam, rtree_dat, rtree_idx, rtree_prefix): # pylint: disable=no-self-use
        """
        This function take the sam file, output of bwa
        and the Rtree_files, and output a baitmap file
        Parameters:
        -----------
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
        copy(rtree_idx, rtree_prefix+".idx")
        copy(rtree_dat, rtree_prefix+".dat")

        idx = index.Rtree(rtree_prefix)

        baitmap = []

        with open(sam_file, "r") as file_in:
            for line in file_in:
                if line[0] != "@":
                    line = line.rstrip().split("\t")

                    try:
                        crm = int(line[2][3:])
                    except IOError:
                        continue

                    srt_pos = int(line[3])
                    end_pos = srt_pos + int(len(line[9]))

                    id_object = idx.intersection(
                        (srt_pos, crm, end_pos, crm),
                        objects=True)

                    hits = [[i.id, i.bbox] for i in id_object]

                    if len(hits) > 1:
                        logger.warning("probe map to two RE fragmnets, " +
                                       " ".join(line)+" start pos"+ str(srt_pos) +
                                       " end pos"+ str(end_pos))

                    elif not hits:
                        logger.warn("Sequence does not"+
                                    "match with any RE fragment, "+
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

        return baitmap


    @task(returns=bool, baitmap_list=IN,
          out_baitmap=FILE_OUT)
    def create_baitmap(self, baitmap_list, out_baitmap): # pylint: disable=no-self-use
        """
        This function takes a list with RE fragments that
        correspond to baits and print it to a file

        Parameters:
        -----------
        baitmap_list: list
            lsit with all the RE fragments corresponding
            to baits
        out: str
            entire pat and name of the .baitmap file
        """
        #print(out_baitmap)
        with open(out_baitmap, "a") as file_out:
            for frag_coord in baitmap_list:
                print("{}\t{}\t{}\t{}\t{}".format(
                    frag_coord[0],
                    frag_coord[1],
                    frag_coord[2],
                    frag_coord[3],
                    "NaN"), file=file_out)

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

        input_bwa = {
            "genome": input_files["genome_fa"],
            "index": input_files["genome_idx"],
            "loc": input_files["probes_fa"]
        }

        output_bwa = {
            "output": output_files["out_bam"]
        }
        metadata_bwa = {
            "genome": Metadata(
                "Assembly", "fasta", input_files["genome_fa"], None,
                {"assembly": "test"}),
            "index": Metadata(
                "index_bwa", "", input_files["genome_fa"],
                {
                    "assembly": "test",
                    "tool": "bwa_indexer"
                }
            ),
            "loc": Metadata(
                "data_chip_seq", "fastq", output_files["out_bam"], None,
                {"assembly": "test"}
            )
        }

        bwa_t = bwaAlignerMEMTool()
        bwa_files, bwa_meta = bwa_t.run(input_bwa, metadata_bwa, output_bwa)

        bwa_meta = compss_wait_on(bwa_meta)

        if "".join(input_files["Rtree_file_dat"].split(".")[:-1]) != \
           "".join(input_files["Rtree_file_idx"].split(".")[:-1]):
            logger.fatal("Rtree_file_dat and Rtree_file_idx"
                         "should have the same prefix name")
            return False

        prefix_rtree = "".join(input_files["Rtree_file_idx"].split(".")[-1])

        baitmap_list = self.sam_to_baitmap(
            output_files["bait_sam"],
            output_files["out_bam"],
            input_files["Rtree_file_dat"],
            input_files["Rtree_file_idx"],
            prefix_rtree)

        results = self.create_baitmap(
            baitmap_list,
            output_files["out_baitmap"])

        results = compss_wait_on(results)

        output_metadata = {
            "out_baitmap": Metadata(
                data_type="RE sites with baits",
                file_type=".baitmap",
                file_path=output_files["out_baitmap"],
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["probes_fa"].file_path,
                    input_metadata["Rtree_file_dat"].file_path,
                    input_metadata["Rtree_file_dat"].file_path,
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "RE" : input_metadata["Rtree_file_idx"].meta_data,
                    "tool": "makeBaitmap",
                }
            ),
            "bait_sam": Metadata(
                data_type="TXT",
                file_type=".sam",
                file_path=output_files["bait_sam"],
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["probes_fa"].file_path,
                    input_metadata["Rtree_file_idx"].file_path,
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "RE" : input_metadata["Rtree_file_idx"].meta_data,
                    "tool": "makeBaitmap",
                }
            ),
            "out_bam": Metadata(
                data_type=".bam",
                file_type=".bam",
                file_path=output_files["out_bam"],
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["probes_fa"].file_path,
                    input_metadata["Rtree_file_idx"].file_path,
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "RE" : input_metadata["Rtree_file_idx"].meta_data,
                    "tool": "makeBaitmap",
                }
            ),

        }

        return output_files, output_metadata
