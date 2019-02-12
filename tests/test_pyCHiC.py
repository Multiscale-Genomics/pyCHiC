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
import pytest # pylint: disable=unused-import
import pandas as pd
from basic_modules.metadata import Metadata
from CHiC.tool.run_pyCHiC import run_pyCHiC

def test_pychic():
    """
    Function to test the read_samples function from pyCHiC
    """
    path = "tests/data/test_run_chicago/data_chicago/"

    input_files = {
        "RMAP" : path +"h19_chr20and21.rmap",
        "BAITMAP" : path +"h19_chr20and21.baitmap",
        "nbpb" : path +"h19_chr20and21.nbpb",
        "npb" : path +"h19_chr20and21.npb",
        "poe" : path +"h19_chr20and21.poe",
        "chinput" :
            path + "GM_rep1.chinput"
            #path + "GM_rep2.chinput",
            #path + "GM_rep3.chinput"
    }

    configuration = {
        "pychic_features_plot" : "DEFB125,DEFB126,DEFB128,DEFB129,DEFB132,TRIB3",
        "pychic_binsize" : 20000,
        "execution" : ".",
        "pychic_cpu" : 3,
        "pychic_cutoff" : 5,
        "pychic_export_format" : ["washU_text"],
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

    metadata = {
        "chinput" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606)
        }

    output_files = {
        "washU_text" : "out_test_washU_text.txt",
        "pdf_examples" : "out_test_examples.pdf",
        "params_out" : "parameters_out.txt",
    }

    pychic_obj = run_pyCHiC(configuration)
    pychic_obj.run(input_files, metadata, output_files)

    output_loc = "out_test_washU_text.txt"
    current_out = pd.read_csv(output_loc, sep="\t")

    output_washu = pd.read_csv("tests/data/output_pychic.txt",
                               sep="\t")

    pd.testing.assert_frame_equal(output_washu,
                                  current_out,
                                  check_less_precise=True,
                                  check_exact=False)
    #remove the files
    import os
    os.remove(configuration["execution"]+"/"+"out_test_washU_text.txt")
    os.remove(configuration["execution"]+"/"+"parameters_out.txt")
    os.remove(configuration["execution"]+"/"+"out_test_examples.pdf")
