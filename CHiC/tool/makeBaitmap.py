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

    @task(returns=bool, genome_index=FILE_IN, probes_fa=FILE_IN,
          out_sam=FILE_OUT, amb=FILE_IN, ann=FILE_IN, bwt=FILE_IN,
          pac=FILE_IN, sa=FILE_IN)
    def bwa_for_probes(self, amb, ann, bwt, pac, sa, genome_index, probes_fa, out_sam):
        """
        This function run bwa using an index genome and a probes file
        in fasta format. bwa is used as single end and with high
        gap penalty and missmacht score

        Parameters:
        -----------
        genome_index: str
            path to the reference genome with indexed files
            in the same folder. This genome should
            be the same used to generate the .rmap file
        probes_fa: str
            path to probes files in fasta format,
            every sequences representing one
        out: str
            Name of the output file
        Return:
        ------
         bool
        """

        args = " ".join(["bwa", "mem", "-O", "100", "-B", "20",
                         genome_index, probes_fa, ">", out_sam])

        logger.info("bwa_for_probes CMD: " + args)

        process = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        if os.path.getsize(out_sam) > 0:
            return True

        logger.fatal("bwa failed to generate sam file")
        logger.fatal("bwa stdout" + proc_out)
        logger.fatal("bwa stderr" + proc_err)
        return False

    @task(returns=bool, genome_index=FILE_IN, probes_fa=FILE_IN,
           genome_fa=FILE_IN, out_sam=FILE_OUT, out_bam=FILE_OUT,
           out_sam=FILE_OUT)
    def bwa_for_probes2(self, genome_fa, genome_index, probes_fa, out_bam, out_sam):
        """
        This function run bwa using an index genome and a probes file
        in fasta format. bwa is used as single end and with high
        gap penalty and missmacht score

        Parameters:
        -----------
        genome_index: str
            path to the reference genome with indexed files
            in the same folder. This genome should
            be the same used to generate the .rmap file
        probes_fa: str
            path to probes files in fasta format,
            every sequences representing one
        out: str
            Name of the output file
        Return:
        ------
         bool
        """

        out_bam = out_sam.split(".")[0]+".bam"

        input_mem = {
            "genome": genome_fa,
            "index": genome_index,
            "loc": probes_fa,
            "bam_loc": out_bam
            }

        output_mem = {
            "output": out_bam,
        }

        metadata_mem = {

            "genome": Metadata("Assembly", "fasta", genome_fa, None, {"assembly": "test"}),

            "index": Metadata("index_bwa", "", [genome_fa], {
                "assembly": "test",
                "tool": "bwa_indexer"
                }
                             ),

            "loc": Metadata("probes", "fastq", probes_fa, None, {"assembly": "test"})
        }

        bwa_aligner = bwaAlignerMEMTool()

        bwa_aligner.run(input_mem, metadata_mem, output_mem)


        tmp_bam = "/".join(probes_fa.split("/")[:-1]) + "/tmp/" + probes_fa.split("/")[-1]+".bam"

        args = ["samtools", "view", "-h", "-o", out_sam, tmp_bam]

        logger.info("samtools args: " + ' '.join(args))

        try:
            with open(out_sam, "w") as f_out:
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

    @task(returns = str, sam_file=FILE_IN, Rtree_files=FILE_IN)
    def sam_to_baitmap(self, sam_file, Rtree_files):
        """
        This function take the sam file, output of bwa
        and the Rtree_files, and output a baitmap file
        Parameters:
        -----------
        sam_file : str
            path to output file from bwa_for_probes
        rmap: str
            complete path to .rmap file
        """
        idx = index.Rtree(Rtree_files)

        baitmap = []

        with open(sam_file, "r") as file_in:
            for line in file_in:
                if line[0] != "@":
                    line = line.rstrip().split("\t")

                    try:
                        crm = int(line[2][3:])
                    except:
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
        return baitmap

    def create_baitmap(self, baitmap_list, out_baitmap):
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
        print(out_baitmap)
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

    def run(self, input_files, metadata, output_files):
        """
        The main function to produce a .baitmap file, starting from rtree files,
        indexed genome and probes.

        Parameters
        ----------
        input_files : Dict
            genome_fa
            genome_idx
            probes_fa
            Rtree_files


        metadata : dict
        Returns
        -------
        output_files : list
            List of locations for the output files.
        output_metadata : list
            List of matching metadata dict objects
        """
        self.bwa_for_probes(
            input_files[".amb"],
            input_files[".ann"],
            input_files[".bwt"],
            input_files[".pac"],
            input_files[".sa"],
            input_files["genome_fa"],
            input_files["genome_idx"],
            input_files["probes_fa"],
            #output_files["out_bam"],
            output_files["out_sam"]
            )

        baitmap_list = self.sam_to_baitmap(
            output_files["out_sam"],
            input_files["Rtree_files"],)

        results = self.create_baitmap(
            baitmap_list,
            output_files["out_baitmap"])

        output_metadata = {
            "baitmap": Metadata(
                data_type="RE sites with baits",
                file_type=".baitmap",
                file_path=output_files["out_baitmap"],
                sources=[
                    metadata["genome_digest"].file_path,
                    metadata["probes"].file_path,
                    metadata["Rtree_files"].file_path,
                ],
                taxon_id=metadata["genome_digest"].taxon_id,
                meta_data={
                    "RE" : metadata["Rtree_files"].meta_data,
                    "tool": "makeBaitmap",
                }
            )
        }


        return results, output_metadata
