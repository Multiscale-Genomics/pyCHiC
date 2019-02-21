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
import sys
import os

from utils import logger
from CHiC.tool.pyCHiC import pyCHiC

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN, OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_delete_file
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN, OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on, compss_delete_file  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------

class run_pyCHiC(Tool):  # pylint: disable=invalid-name
    """
    Tool for runnning pyCHiC
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("starting run_pyCHiC")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)


    @task(returns=bool, configuration=IN, RMAP=FILE_IN, BAITMAP=FILE_IN,
          npb=FILE_IN, nbpb=FILE_IN, poe=FILE_IN, chinput=FILE_IN)
    def pychic_runner(self, configuration, RMAP, BAITMAP, npb, nbpb, poe, chinput,
                      washU_text, pdf_examples, params_out):
        """
        This is the function with pyCOMPSs that is going to run pyCHiC

        Parameters
        ----------

        Returns
        -------
        """
        input_files = {"RMAP": RMAP,
                       "BAITMAP": BAITMAP,
                       "npb": npb,
                       "nbpb": nbpb,
                       "poe": poe,
                       "chinput": chinput
                       }

        metadata = {}

        output_files = {"washU_text": washU_text,
                        "pdf_examples": pdf_examples,
                        "params_out": params_out}

        pychic_handler = pyCHiC(configuration)
        pychic_handler.run(input_files, metadata, output_files)

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run the test_writer tool

        Parameters
        ----------
        input_files : dict
            List of input files - In this case there are no input files required
        input_metadata: dict
            Matching metadata for each of the files, plus any additional data
        output_files : dict
            List of the output files that are to be generated

        Returns
        -------
        output_files : dict
            List of files with a single entry.
        output_metadata : dict
            List of matching metadata for the returned files
        """

        rtree_file_dat = "tests/data/test_rmap/rtree_file.dat"
        rtree_file_idx = "tests/data/test_rmap/rtree_file.idx"
        chr_handler = "tests/data/test_baitmap/chr_handler.txt"
        rmap = "tests/data/test_run_chicago/test.rmap"
        baitmap = "tests/data/test_run_chicago/test.baitmap"
        bait_sam = "tests/data/test_baitmap/baits.sam"
        nbpb = os.path.abspath("tests/data/test_run_chicago/test.nbpb")
        npb = os.path.abspath("tests/data/test_run_chicago/test.npb")
        poe = os.path.abspath("tests/data/test_run_chicago/test.poe")
        out_bam = "tests/data/test_baitmap/baits.bam"
        sorted_bam = self.configuration["execution"] + "/" + "sorted_bam"

        if "RMAP" not in input_files:
            input_files["RMAP"] = rmap
        if "BAITMAP" not in input_files:
            input_files["BAITMAP"] = baitmap
        if "nbpb" not in input_files:
            input_files["nbpb"] = nbpb
        if "npb" not in input_files:
            input_files["npb"] = npb
        if "poe" not in input_files:
            input_files["poe"] = poe
        if "chinput" not in input_files:
            input_files["chinput"] = output_files["chinput"]

        if "pychic_binsize" not in self.configuration:
            self.configuration["pychic_binsize"] = \
                int(self.configuration["makeDesignFiles_binsize"])
        else:
            self.configuration["pychic_binsize"] = int(self.configuration["pychic_binsize"])

        if "pychic_minFragLen" not in self.configuration:
            self.configuration["pychic_minFragLen"] = \
                int(self.configuration["makeDesignFiles_minFragLen"])
        else:
            self.configuration["pychic_minFragLen"] = int(self.configuration["pychic_minFragLen"])

        if "pychic_maxFragLen" not in self.configuration:
            self.configuration["pychic_maxFragLen"] = \
                int(self.configuration["makeDesignFiles_maxFragLen"])
        else:
            self.configuration["pychic_maxFragLen"] = int(self.configuration["pychic_maxFragLen"])

        if "pychic_maxLBrownEst" not in self.configuration:
            self.configuration["pychic_maxLBrownEst"] = \
                float(self.configuration["makeDesignFiles_maxLBrownEst"])
        else:
            self.configuration["pychic_maxLBrownEst"] = \
                float(self.configuration["pychic_maxLBrownEst"])

        self.configuration["pychic_removeAdjacent"] = True
        self.configuration["pychic_adjBait2bait"] = True

        if "pychic_bam" not in self.configuration:
            self.configuration["pychic_bam"] = sorted_bam

        pychic_handler = pyCHiC(self.configuration)
        pychic_handler.run(input_files, input_metadata, output_files)

        #results = self.pychic_runner(self.configuration,
        #                             input_files["RMAP"],
        #                             input_files["BAITMAP"],
        #                             input_files["npb"],
        #                             input_files["nbpb"],
        #                             input_files["poe"],
        #                             input_files["chinput"],
        #                             output_files["washU_text"],
        #                             output_files["pdf_examples"],
        #                                                output_files["params_out"])

        if "genome_name" in self.configuration:
            files_dir = os.listdir(self.configuration["execution"])
            for file_ in files_dir:
                if file_.startswith("Digest_"+self.configuration["genome_name"]):
                    os.remove(file_)

        compss_delete_file(rtree_file_idx)
        compss_delete_file(rtree_file_dat)
        compss_delete_file(chr_handler)
        compss_delete_file(rmap)
        compss_delete_file(baitmap)
        compss_delete_file(bait_sam)
        compss_delete_file(npb)
        compss_delete_file(nbpb)
        compss_delete_file(poe)
        compss_delete_file(out_bam)
        compss_delete_file(sorted_bam)

        if "chinput" not in input_metadata:
            input_metadata["chinput"] = input_metadata["genome_fa"]

        output_metadata = {
            "washU_text" : Metadata(
                data_type="data_chic",
                file_type="TXT",
                file_path=output_files["washU_text"],
                sources=[

                ],
                taxon_id=input_metadata["chinput"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago"
                }
            ),

            "pdf_examples" : Metadata(
                data_type="data_chic",
                file_type="PDF",
                file_path=output_files["pdf_examples"],
                sources=[
                ],
                taxon_id=input_metadata["chinput"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago"
                }
            ),

            "params_out" : Metadata(
                data_type="data_chic",
                file_type="TXT",
                file_path=output_files["params_out"],
                sources=[
                ],
                taxon_id=input_metadata["chinput"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago"
                }
            )
        }

        return output_files, output_metadata

"""
if __name__ == "__main__":

    path = "../../tests/data/test_run_chicago/data_chicago/"

    input_files = {
        "RMAP" : path +"h19_chr20and21.rmap",
        "BAITMAP" : path +"h19_chr20and21.baitmap",
        "nbpb" : path +"h19_chr20and21.nbpb",
        "npb" : path +"h19_chr20and21.npb",
        "poe" : path +"h19_chr20and21.poe",
        "chinput" : path + "GM_rep1.chinput"
    }

    metadata = {
        "chinput" : Metadata(
        "data_chicago", "chinput", [], None, None, 9606)
    }

    output_files = {
        "washU_text" : "out_test_washU_text.txt",
        "pdf_examples" : "out_test_examples.pdf",
        "params_out" : "parameters.txt"
        }

    configuration = {
        "pychic_features_plot" : None,
        "pychic_binsize" : 20000,
        "execution" : ".",
        "pychic_cpu" : 3,
        "pychic_cutoff" : 5,
        "pychic_export_format" : ["washU_text", "seqMonk", "interBed"],
        "pychic_order" : "score",
        "pychic_maxLBrownEst" : 1500000.0,
        "pychic_minFragLen" : 150, # minimun OE fragment lenght in bps
        "pychic_maxFragLen" : 40000, # maximun OE fragment lenght in bps
        "pychic_minNPerBait" : 250, # total number of interactions per bait
        "pychic_removeAdjacent" : "True",
        "pychic_adjBait2bait" : "True",
        "pychic_tlb_filterTopPercent" : 0.01,
        "pychic_tlb_minProxOEPerBin" : 150,
        "pychic_tlb_minProxB2BPerBin" : 15,
        "pychic_techNoise_minBaitsPerBin" : 150,
        "pychic_brownianNoise_samples" : 1,
        "pychic_brownianNoise_subset" : 500,
        "pychic_brownianNoise_seed" : 3,
        "pychic_weightAlpha" : 34.1157346557331,
        "pychic_weightBeta" : -2.58688050486759,
        "pychic_weightGamma" : -17.1347845819659,
        "pychic_weightDelta" : -7.07609245521541,
        #"pychic_Rda" : "False",
        #"pychic_output_dir" : path +"output_pyCHiC",

    }

    pyCHiC_obj = run_pyCHiC(configuration)
    pyCHiC_obj.run(input_files, metadata, output_files)
"""