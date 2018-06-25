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
from shutil import copy2

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
from tool.aligner_utils import alignerUtils

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


    @task(returns=bool, genome_file_name=IN, genome_idx=FILE_IN,
          amb_file=FILE_OUT, ann_file=FILE_OUT, bwt_file=FILE_OUT,
          pac_file=FILE_OUT, sa_file=FILE_OUT)
    def untar_index(  # pylint: disable=too-many-locals,too-many-arguments
            self, genome_file_name, genome_idx,
            amb_file, ann_file, bwt_file, pac_file, sa_file):
        """
        Extracts the BWA index files from the genome index tar file.
        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        genome_idx : str
            Location of the BWA index file
        amb_file : str
            Location of the amb index file
        ann_file : str
            Location of the ann index file
        bwt_file : str
            Location of the bwt index file
        pac_file : str
            Location of the pac index file
        sa_file : str
            Location of the sa index file
        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            return True

        gfl = genome_file_name.split("/")
        au_handle = alignerUtils()
        au_handle.bwa_untar_index(
            gfl[-1], genome_idx, amb_file, ann_file, bwt_file, pac_file, sa_file)

        return True


    @task(returns=bool, amb=FILE_IN, ann=FILE_IN, bwt=FILE_IN,
          pac=FILE_IN, sa=FILE_IN, genome_index=FILE_IN, probes_fa=FILE_IN,
          bait_sam=FILE_OUT, )
    def bwa_for_probes(self, amb, ann, bwt, pac, sa, genome_index, probes_fa, bait_sam):
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
                         genome_index, probes_fa, ">", bait_sam])

        logger.info("bwa_for_probes CMD: " + args)

        process = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        if os.path.getsize(bait_sam) > 0:
            return True

        logger.fatal("bwa failed to generate sam file")
        logger.fatal("bwa stdout" + proc_out)
        logger.fatal("bwa stderr" + proc_err)
        return False

    @task(returns=bool, genome_index=FILE_IN, probes_fa=FILE_IN,
          genome_fa=FILE_IN, bait_sam=FILE_OUT, out_bam=FILE_OUT,
          amb=FILE_IN, ann=FILE_IN, bwt=FILE_IN, pac=FILE_IN,
          sa=FILE_IN)
    def bwa_for_probes2(self, amb, ann, bwt, pac, sa, genome_fa,
                        probes_fa, out_bam, bait_sam):
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

        out_bam = bait_sam.split(".")[0]+".bam"

        input_mem = {
            "genome": genome_fa,
            "index": genome_fa,
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

        compss_wait_on(bwa_aligner)

        tmp_bam = "/".join(probes_fa.split("/")[:-1]) + "/tmp/" + probes_fa.split("/")[-1]+".bam"

        args = ["samtools", "view", "-h", "-o", bait_sam, tmp_bam]

        logger.info("samtools args: " + ' '.join(args))


        try:
            with open(bait_sam, "w") as f_out:
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

    @task(returns = str, sam_file=FILE_IN, rtree_dat=FILE_IN, rtree_idx=FILE_IN,
          rtree_prefix=IN)
    def sam_to_baitmap(self, sam_file, rtree_dat, rtree_idx, rtree_prefix):
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
        copy2(rtree_idx, rtree_prefix+".idx")
        copy2(rtree_dat, rtree_prefix+".dat")

        idx = index.Rtree(rtree_prefix)

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

        os.remove(rtree_prefix+".dat")
        os.remove(rtree_prefix+".idx")

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

        index_files = {
            "amb": input_files["genome_fa"] + ".amb",
            "ann": input_files["genome_fa"] + ".ann",
            "bwt": input_files["genome_fa"] + ".bwt",
            "pac": input_files["genome_fa"] + ".pac",
            "sa": input_files["genome_fa"] + ".sa"
        }

        self.untar_index(input_files["genome_fa"],
                         input_files["genome_idx"],
                         index_files["amb"],
                         index_files["ann"],
                         index_files["bwt"],
                         index_files["pac"],
                         index_files["sa"],
                         )

        self.bwa_for_probes2(
            index_files["amb"],
            index_files["ann"],
            index_files["bwt"],
            index_files["pac"],
            index_files["sa"],
            input_files["genome_fa"],
            input_files["probes_fa"],
            output_files["out_bam"],
            output_files["bait_sam"]
            )

        if "".join(input_files["Rtree_file_dat"].split(".")[:-1]) != \
           "".join(input_files["Rtree_file_idx"].split(".")[:-1]):
            logger.fatal("Rtree_file_dat and Rtree_file_idx"
                         "should have the same prefix name")
            return False

        prefix_rtree = "".join(input_files["Rtree_file_idx"].split(".")[-1])

        baitmap_list = self.sam_to_baitmap(
            output_files["bait_sam"],
            input_files["Rtree_file_dat"],
            input_files["Rtree_file_idx"],
            prefix_rtree,)

        results = self.create_baitmap(
            baitmap_list,
            output_files["out_baitmap"])

        output_metadata = {
            "baitmap": Metadata(
                data_type="RE sites with baits",
                file_type=".baitmap",
                file_path=output_files["out_baitmap"],
                sources=[
                    metadata["genome_fa"].file_path,
                    metadata["probes"].file_path,
                    metadata["Rtree_file_dat"].file_path,
                    metadata["Rtree_file_dat"].file_path,
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "RE" : metadata["Rtree_file_idx"].meta_data,
                    "tool": "makeBaitmap",
                }
            ),
            "bait_sam": Metadata(
                data_type="TXT",
                file_type=".sam",
                file_path=output_files["bait_sam"],
                sources=[
                    metadata["genome_fa"].file_path,
                    metadata["probes"].file_path,
                    metadata["Rtree_file_idx"].file_path,
                ],
                taxon_id=metadata["genome_fa"].taxon_id,
                meta_data={
                    "RE" : metadata["Rtree_file_idx"].meta_data,
                    "tool": "makeBaitmap",
                }
            )
        }

        return results, output_metadata

if __name__ == "__main__":

    path = "../../tests/data/"

    configuration = {
    }

    input_files = {
        "genome_idx" : path + "test_baitmap/chr21_hg19.fa.tar.gz",
        "probes_fa" : path + "test_baitmap/baits.fa",
        "Rtree_file_dat" : path + "test_rmap/rtree_file.dat",
        "Rtree_file_idx" : path + "test_rmap/rtree_file.idx",
        "genome_fa" : path+ "test_baitmap/chr21_hg19.fa"
    }

    output_files = {
        "bait_sam" :  path + "test_baitmap/baits.sam",
        "out_bam" : path +  "tests/baits.bam",
        "out_baitmap" : path + "test_run_chicago/test.baitmap"
    }

    metadata = {
        "genome_fa" : Metadata(
            "hg38", "fasta", path + "test_rmap/chr21_hg19.fa",
            None, "HindIII", 9606),

        "probes" : Metadata(
            "C-HiC probes", "fasta", path + "test_baitmap/baits.fa",
            None, None, 9606),

        "Rtree_file_dat" : Metadata(
            "Rtree files", [".dat", ".idx"], path + "test_rmap/rtree_file",
            {"genome" : path + "test_rmap/chr21_hg19.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606
            ),

        "Rtree_file_idx" : Metadata(

            "Rtree files", [".dat", ".idx"], path + "test_rmap/rtree_file",
            {"genome" : path + "test_rmap/chr21_hg19.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606
            )
    }

    baitmap_handler = makeBaitmapTool(configuration)
    baitmap_handler.run(input_files, metadata, output_files)