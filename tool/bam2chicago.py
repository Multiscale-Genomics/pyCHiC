
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

import argparse
import sys
import numpy as np
import pandas as pd
import subprocess


class bam2chicago(Tool):

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


    def check_files(self, RMAP, BAITMAP):
        """
        This function check that the rmap and baitmap input fiels are in the correct format

        PARAMETERS:
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
        except:
            logger.fatal("rmap rows contain"+
                         "different number of columns")

        if rmap.shape[1] != 4:
            logger.fatal("rmap file does not have 4 columns")

        rmap_IDs = rmap.iloc[:,3]
        set_rmap_IDs = set(rmap_IDs)

        if len(rmap_IDs) != len(set_rmap_IDs):
            logger.fatal("Error! Duplicated fragment IDs found in rmap file")

        print("Checking baitmap file")
        try:
            baitmap = pd.read_table(BAITMAP, header=None)
        except:
            logger.fatal("baitmap rows contain"+
                         "different number of columns")

        if baitmap.shape[1] != 4:
            logger.fatal("baitmap does not contain 4 columns")

        baitmapIDs = baitmap.iloc[:, 3]
        set_baitmapIDs = set(baitmapIDs)

        if len(baitmapIDs) != len(set_baitmapIDs):
            logger.fatal("Error! Duplicated fragment IDs found in baitmap")

        baits_no_rmap = set_baitmapIDs - set_rmap_IDs

        if not baits_no_rmap:
            print("rmap and baitmap files checked successfully =)")
        else:
            logger.fatal("Error! IDs Baitmap entry not found in rmap: " + baits_no_rmap)

        return True

    def intersect_bam_bait(self, BAM, BAITMAP, RMAP, out_dir):
        """
        THis function will intersect the pair-end reads from the
        bam file with the baits
        PARAMETERS:
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
                       BAITMAP, "-f", "0.6",">", out_dir+"/mappedToBaits.bedpe"]

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

        logger.info("Flipping all reads that overlap with the bait on to the right-hand side...")

        with open(out_dir+"/mappedToBaits.bedpe", "r") as FILE_IN:
            with open(out_dir+"/mappedToBaits_baitOnRight.bedpe", "w") as FILE_OUT:
                for line in FILE_IN:
                    line_hdl = line.rstrip().split("\t")
                    minRight=min(int(line_hdl[2]), int(line_hdl[12]))
                    maxLeft=max(int(line_hdl[1]), int(line_hdl[11]))

                    if (line_hdl[0] == line_hdl[10]) and (minRight - maxLeft) / (int(line_hdl[2]) - int(line_hdl[1])) >= 0.6:
                        print("\t".join([line_hdl[3], line_hdl[4], line_hdl[5], line_hdl[0], line_hdl[1], line_hdl[2],
                            line_hdl[6], line_hdl[7], line_hdl[9], line_hdl[8], line_hdl[10], line_hdl[11],
                            line_hdl[12], line_hdl[13]]), file=FILE_OUT)
                    else:
                        print("\t".join(line_hdl), file=FILE_OUT)

        logger.info("Intersecting with bait fragments again to produce a list of bait-to-bait"+
                    "interactions that can be used separately; note they will also be retained"+
                    "in the main output...")

        with open(out_dir+"/bait2bait.bedpe", "w") as b2b:
            print("## samplename="+out_dir+" bamname="+BAM+ " baitmap="+BAITMAP+ " rmap="+RMAP, file=b2b)

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

        logger.info("bedtools inteersect parameters: "
         + " ".join(bashCommand))

        process = subprocess.Popen(" ".join(basCommand), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr= subprocess.PIPE)
        process.wait()
        output, error = process.communicate()

        return True

    def filter_reads(self, out_dir, mappedToBaitsBoRAndRFrag):
        """
        This function filter reads that do not pass the minimun required overhand.
        PARAMTERS:
        ----------
        out_dir: str
            path to the ouput directoruy
        mappedToBaitsBoRAndRFrag: str
            name of the output generated by bedtools intersection with
            restriction enzyme.
        RETURNS:
        -------
        bool
        """
        mapped_reads = pd.read_csv(mappedToBaitsBoRAndRFrag, sep= '\t',
                                   header=None)

        less06 = mapped_reads.loc[(mapped_reads.iloc[:,15] == -1) & \
                                  (mapped_reads.iloc[:,16] == -1)]

        less06.to_csv(out_dir+"/mappedToBaitsBoRAndRFrag_fless06.bedpe", \
                      sep="\t", index=False, header=None)

        less06_rows = less06.shape[0]

        more06 = mapped_reads.loc[(mapped_reads.iloc[:,15] != -1) & \
                                  (mapped_reads.iloc[:,16] != -1)]

        more06.to_csv(out_dir+"/mappedToBaitsBoRAndRFrag_more062.bedpe", \
                      sep="\t", index=False, header=None)

        more06_rows = more06.shape[0]

        logger.info("Filtered out %f reads with <60%% overlap with a \
                    single digestion fragment" % \
                    (less06_rows/(less06_rows+more06_rows)))

        return True

    def calculate_distances(self, out_dir, mappedToBaitsBoRAndRFrag_fmore06):
        """
        adasda
        PARAMETERS:
        -----------
        mappedToBaitsBoRAndRFrag_fmore06: str
            path to the file

        """
        logger.info("Adding frag length and signed distance from bait; "\
                    "removing self-ligation fragments (if any; not expected" \
                    "with HiCUP input)...")

        with open(mappedToBaitsBoRAndRFrag_fmore06, "r") as FILE_IN:
            with open(out_dir+"/mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe", "w") as FILE_OUT:
                for line in FILE_IN:
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
                            print("\t".join(line_hdl), file = FILE_OUT)


        if os.path.getsize(out_dir+"//mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe") > 0:
            pass
        else:
            logger.fatal("calculate_distances function failed")

    def write_output(self, out_dir, RMAP, mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen):
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

        with open(mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen, "r") as FILE_IN:
            with open(out_dir+"/"+out_dir+"unsorted.chinput", "w") as FILE_OUT:
                print("##samplename'\t'bamname='\t'baitmapfile='\t'digestfile='\t'")
                for line in FILE_IN:
                    line_hdl = line.rstrip().split("\t")
                    if line_hdl[13]+" "+line_hdl[17] not in baits_dict:
                        baits_dict[line_hdl[13]+" "+line_hdl[17]] = [1, line_hdl[19], line_hdl[20]]
                    else:
                        baits_dict[line_hdl[13]+" "+line_hdl[17]] = [baits_dict[line_hdl[13]+" "+line_hdl[17]][0] + 1,
                                                                     baits_dict[line_hdl[13]+" "+line_hdl[17]][1],
                                                                     baits_dict[line_hdl[13]+" "+line_hdl[17]][2]]
                for item in baits_dict:
                    print(item.split(" ")[0]+"\t"+item.split(" ")[1]+
                          "\t"+"\t".join(list(map(lambda x : str(x),
                          baits_dict[item]))), file= FILE_OUT)


        bashCommand = ["cat", out_dir+"/"+out_dir+"unsorted.chinput", "|" ,\
                      "sort", "-k1,1", "-k2,2n", "-T",RMAP,">>",
                      out_dir+"/"+out_dir+".chinput"]

        logger.info("Sorting output: "
         + " ".join(bashCommand))

        process = subprocess.Popen(" ".join(bashCommand), shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr= subprocess.PIPE)



if __name__ == "__main__":

    test = bam2chicago()
    """
    test.intersect_bam_bait("/Users/pacera/developing/C-HiC/genomes/out_hicup/SRR3535023_1_2.hicup.bam",
        "/Users/pacera/developing/C-HiC/tests/data/hg19TestDesign/h19_chr20and21.baitmap_4col_chr.txt",
        "/Users/pacera/developing/C-HiC/tests/data/hg19TestDesign/h19_chr20and21_chr.rmap",
        "output_test")
    """

    #test.filter_reads("output_test", "/Users/pacera/developing/C-HiC/tool/output_test/mappedToBaitsBoRAndRFrag.bedpe")

    test.write_output("output_test", "/Users/pacera/developing/C-HiC/tests/data/hg19TestDesign/h19_chr20and21.rmap"  ,"/Users/pacera/developing/C-HiC/tool/output_test/mappedToBaitsBoRAndRFrag_fmore06_withDistSignLen.bedpe")










































