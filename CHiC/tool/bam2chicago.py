
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
from __future__ import division
from __future__ import print_function

import sys
import os
import subprocess
import pandas as pd

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



class bam2chicago(Tool):

    """
    Tool for digest the genome with one RE. Wrapper of hicup_digester
    """

    def __init__(self, configuration=None):
        """
        initialising the class

        Parameters
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


    def check_files(self, RMAP, BAITMAP, BAM):
        """
        This function check that the rmap and baitmap input fiels are in the correct format

        Parameters
        -----------
        RMAP: str
            path to rmap file, 4 columns
        BAITMAP: str
            path to baitmap file, 4 columns

        RETURN:
        -------
        Bool

        """
        print("Checking rmap file")
        try:
            rmap = pd.read_table(RMAP, header=None)
        except IOError:
            logger.fatal("rmap rows contain"+
                         "different number of columns")
            return False

        if rmap.shape[1] != 4:
            logger.fatal("rmap file does not have 4 columns")

        #retrieve rmap IDs and chromosome format
        rmap_chr = rmap.iloc[1, 1]
        rmap_IDs = rmap.iloc[:, 3]
        set_rmap_IDs = set(rmap_IDs)

        if len(rmap_IDs) != len(set_rmap_IDs):
            logger.fatal("Error! Duplicated fragment IDs found in rmap file")

        print("Checking baitmap file")
        try:
            baitmap = pd.read_table(BAITMAP, header=None)
        except IOError:
            logger.fatal("baitmap rows contain"+
                         "different number of columns")
            return False

        if baitmap.shape[1] != 4:
            print(baitmap.shape)
            logger.fatal("baitmap does not contain 4 columns")

        #retrieve baitmap IDs
        baitmap_chr = baitmap.iloc[1, 1]
        baitmapIDs = baitmap.iloc[:, 3]
        set_baitmapIDs = set(baitmapIDs)

        if len(baitmapIDs) != len(set_baitmapIDs):
            logger.fatal("Error! Duplicated fragment IDs found in baitmap")

        baits_no_rmap = set_baitmapIDs - set_rmap_IDs

        if baits_no_rmap:
            logger.fatal("Error! IDs Baitmap entry not found in rmap: " + baits_no_rmap)
            return False
        else:
            logger.info("rmap and baitmap files checked successfully =)")

        return True

    def intersect_bam_bait_rmap(self, BAM, BAITMAP, RMAP, out_dir):
        """
        THis function will intersect the pair-end reads from the
        bam file with the baits

        Parameters
        -----------
        BAM: str
            path to the bam file

        BAITMAP: str
            path to the bait file
        """
        logger.info("Intersecting with bait fragments (using min overhang of 0.6)...")

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        bashCommand = ["bedtools", "pairtobed", "-abam", BAM, "-bedpe", "-b",
                       BAITMAP, "-f", "0.6", ">", out_dir+"/mappedToBaits.bedpe"]

        logger.info("bedtools pairtobed parameters :"+ " ".join(bashCommand))
        process = subprocess.Popen(" ".join(bashCommand), shell=True,
                                   stdout=subprocess.PIPE)
        process.wait()
        output, error = process.communicate()

        if os.path.getsize(out_dir+"/mappedToBaits.bedpe") > 0:
            pass
        else:
            logger.fatal("bedtools failed in this operation" \
                "parameters :"+ " ".join(bashCommand))
            logger.fatal(output)
            logger.fatal(error)
            logger.fatal("Check that the chromosome format is the same" \
                         " in rmap, baitmap and bam files")

        logger.info("Flipping all reads that overlap with the bait on to the right-hand side...")

        with open(out_dir+"/mappedToBaits.bedpe", "r") as file_in:
            with open(out_dir+"/mappedToBaits_baitOnRight.bedpe", "w") as file_out:
                for line in file_in:
                    line_hdl = line.rstrip().split("\t")
                    minRight = min(int(line_hdl[2]), int(line_hdl[12]))
                    maxLeft = max(int(line_hdl[1]), int(line_hdl[11]))

                    if (line_hdl[0] == line_hdl[10]) and \
                       (minRight - maxLeft) / (int(line_hdl[2]) - int(line_hdl[1])) >= 0.6:
                        print("\t".join([line_hdl[3], line_hdl[4], line_hdl[5],
                              line_hdl[0], line_hdl[1], line_hdl[2],
                              line_hdl[6], line_hdl[7], line_hdl[9], line_hdl[8],
                              line_hdl[10], line_hdl[11],
                              line_hdl[12], line_hdl[13]]), file=file_out)
                    else:
                        print("\t".join(line_hdl), file=file_out)

        logger.info("Intersecting with bait fragments again to produce a list of bait-to-bait"+
                    "interactions that can be used separately; note they will also be retained"+
                    "in the main output...")

        with open(out_dir+"/bait2bait.bedpe", "w") as b2b:
            print("## samplename="+out_dir+" bamname="+BAM+ \
                  " baitmap="+BAITMAP+ " rmap="+RMAP, file=b2b)

        bashCommand = ["bedtools", "intersect", "-a", out_dir+"/mappedToBaits_baitOnRight.bedpe",
                       "-wo", "-f", "0.6", "-b", BAITMAP, ">>", out_dir+ "/bait2bait.bedpe"]

        logger.info("bedtools intersect parameters: " + " ".join(bashCommand))
        process = subprocess.Popen(" ".join(bashCommand), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        output, error = process.communicate()

        if os.path.getsize(out_dir+"/bait2bait.bedpe") > 0:
            pass
        else:
            logger.fatal("bedtools failed in this operation" \
                "parameters :"+ " ".join(bashCommand))


        logger.info("Intersecting with restriction"\
            "fragments (using min overhang of 0.6)...")

        basCommand = ["bedtools", "intersect", "-a",
                      out_dir+"/mappedToBaits_baitOnRight.bedpe",
                      "-wao", "-f", "0.6", "-b", RMAP, ">",
                      out_dir+"/mappedToBaitsBoRAndRFrag.bedpe"]

        logger.info("bedtools intersect parameters: "
                    + " ".join(bashCommand))

        process = subprocess.Popen(" ".join(basCommand), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        output, error = process.communicate()

        return out_dir+"/mappedToBaitsBoRAndRFrag.bedpe"

    def filter_reads(self, out_dir, mappedToBaitsBoRAndRFrag):
        """
        This function filter reads that do not pass the minimun required overhand.

        Parameters
        ----------
        out_dir: str
            path to the ouput directoruy
        mappedToBaitsBoRAndRFrag: str
            name of the output generated by bedtools intersection with
            restriction enzyme.

        Returns
        -------
        bool
        """
        mapped_reads = pd.read_csv(mappedToBaitsBoRAndRFrag, sep= '\t',
                                   header=None)

        less06 = mapped_reads.loc[(mapped_reads.iloc[:, 15] == -1) & \
                                  (mapped_reads.iloc[:, 16] == -1)]

        less06.to_csv(out_dir+"/mappedToBaitsBoRAndRFrag_fless06.bedpe", \
                      sep="\t", index=False, header=None)

        less06_rows = less06.shape[0]

        more06 = mapped_reads.loc[(mapped_reads.iloc[:, 15] != -1) & \
                                  (mapped_reads.iloc[:, 16] != -1)]

        more06.to_csv(out_dir+"/mappedToBaitsBoRAndRFrag_fmore06.bedpe", \
                      sep="\t", index=False, header=None)

        more06_rows = more06.shape[0]

        logger.info("Filtered out %f reads with <60%% overlap with a" \
                    "single digestion fragment" % \
                    (less06_rows/(less06_rows+more06_rows)))

        return out_dir+"/mappedToBaitsBoRAndRFrag_fmore06.bedpe"

    def calculate_distances(self, out_dir, mappedToBaitsBoRAndRFrag_fmore06):
        """
        PARAMETERS:
        -----------
        mappedToBaitsBoRAndRFrag_fmore06: str
            path to the file

        """
        logger.info("Adding frag length and signed distance from bait; "\
                    "removing self-ligation fragments (if any; not expected" \
                    "with HiCUP input)...")

        with open(mappedToBaitsBoRAndRFrag_fmore06, "r") as file_in:
            with open(out_dir+"/mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe",
                      "w") as file_out:

                for line in file_in:
                    line_hdl = line.rstrip().split("\t")
                    if line_hdl[10] == line_hdl[14]:
                        midB = int((int(line_hdl[12]) +
                                    int(line_hdl[11]))/2+0.5)

                        midR = int((int(line_hdl[15]) +
                                    int(line_hdl[16]))/2+0.5)
                        distance = midR - midB
                        if distance != 0:
                            line_hdl.append(str(int(line_hdl[16]) - int(line_hdl[15])))
                            line_hdl.append(str(distance))
                            print("\t".join(line_hdl), file=file_out)


        if os.path.getsize(out_dir+"/mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe") > 0:
            pass
        else:
            logger.fatal("calculate_distances function failed")
            return False

        return out_dir+"/mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe"

    def write_output(self, out_dir, RMAP, BAM, BAITMAP, output_file,
                     mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen,):
        """
        This function writes the output in the chinput format
        PARAMETERS:
        ----------
        mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen: str
            path to the file

        RETURN:
        ------
        bool
        """
        logger.info("Pooling read pairs...")

        baits_dict = {}

        with open(mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen, "r") as file_in:
            with open(out_dir+output_file+"unsorted.chinput", "w") as file_out:
                print("##samplename="+out_dir+"'\t'bamname="+BAM+" \
                      '\t'baitmapfile="+BAITMAP+"'\t'rmap="+RMAP+"'\t'", file=file_out)
                for line in file_in:
                    line_hdl = line.rstrip().split("\t")
                    if line_hdl[13]+" "+line_hdl[17] not in baits_dict:
                        baits_dict[line_hdl[13]+" "+line_hdl[17]] = [1, line_hdl[19], line_hdl[20]]
                    else:
                        baits_dict[line_hdl[13]+" "+line_hdl[17]] = [
                            baits_dict[line_hdl[13]+" "+line_hdl[17]][0] + 1,
                            baits_dict[line_hdl[13]+" "+line_hdl[17]][1],
                            baits_dict[line_hdl[13]+" "+line_hdl[17]][2]
                            ]

                for item in baits_dict:
                    print(item.split(" ")[0]+"\t"+item.split(" ")[1]+
                          "\t"+"\t".join([str(i) for i in baits_dict[item]]), file=file_out)

        #sort the printed file based on RMAP file
        bashCommand = ["cat", out_dir+output_file+"unsorted.chinput", "|" ,\
                      "sort", "-k1,1", "-k2,2n", "-T",RMAP,">>",
                      out_dir+output_file+".chinput"]

        logger.info("Sorting output: "
         + " ".join(bashCommand))

        process = subprocess.Popen(" ".join(bashCommand), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr= subprocess.PIPE)
        process.wait()

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        This function runs all the fucntions and produce the file output
        PARAMETERS
        ----------
        input_files: dict
            RMAP: str
                path to rmap file
            BAITMAP: str
                path to the baitmap file
            BAM: str
                path to the bam file

        input_metadata: dict
            input metadata

        output_files:dict
            out_dir: str
                path to the output directory
            output_file: str
                name of the .chinput file
        """

        check = self.check_files(input_files["RMAP"],
                                 input_files["BAITMAP"],
                                 input_files["BAM"])
        if check is True:
            intersected_files = self.intersect_bam_bait_rmap(input_files["BAM"],
                                                             input_files["BAITMAP"],
                                                             input_files["RMAP"],
                                                             output_files["out_dir"])

            filtered_reads = self.filter_reads(output_files["out_dir"],
                                               intersected_files)

            calculated_distances = self.calculate_distances(output_files["out_dir"],
                                                            filtered_reads)

            results = self.write_output(output_files["out_dir"],
                                        input_files["RMAP"],
                                        input_files["BAM"],
                                        input_files["BAITMAP"],
                                        output_files["output_file"],
                                        calculated_distances
                                        )

            output_metadata = {
                ".chinput" : Metadata(
                    data_type="txt",
                    file_type=".chinput",
                    file_path=output_files["out_dir"]+output_files["output_file"],
                    sources=[
                        input_metadata[".rmap"].file_path,
                        input_metadata[".baitmap"].file_path,
                        input_metadata["bam"].file_path
                    ],
                    taxon_id=9606,
                    meta_data={
                        "tool": "Chicago, capture Capture-HiC algorithm"
                    }
                )
            }

        if os.path.getsize(output_files["out_dir"]+output_files["output_file"]+".chinput") > 0:
            pass
        else:
            logger.fatal("write_output function failed to "\
                         "generate output")
        return True

if __name__ == "__main__":

    path = "../test/data/"

    input_files = {
        "RMAP" : path + "test_bam2chicago/chrtest.rmap",
        "BAITMAP" : path +  "test_bam2chicago/chrtest4.baitmap",
        "BAM" : path + "test_bed2bam/outbam_sorted.bam"
    }

    output_files = {
        "out_dir" : path,
        "output_file" :  "test_bam2chicago/out_py"
    }

    input_metadata = {
        ".rmap" : Metadata(
            "data_chicago_input", ".rmap",
            path+"/h19_chr20and21_chr.rmap", None, {}, 9606),
        ".baitmap" : Metadata(
            "data_chicago_input", ".baitmap",
            path+"/h19_chr20and21.baitmap_4col_chr.txt", None, {}, 9606),
        "bam" : Metadata(
            "txt", "bamfile", path + "/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606)
    }

    bam2chicago_handle = bam2chicago()
    bam2chicago_handle.run(input_files, input_metadata, output_files)































