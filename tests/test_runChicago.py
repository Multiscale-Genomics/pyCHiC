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

import os.path
import pytest

from basic_modules.metadata import Metadata
from tool.runChicago import ChicagoTool

def test_chicago():
    """
    Function to testing the R chicago wrapper runChicago.py
    """
    #
    input_files = {
        "chinput_files":
            "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/"+
            "PCHiCdata/inst/extdata/GMchinputFiles/GM_rep1.chinput"
        }

    output_files = {
        "output_dir":
            "/Users/pacera/test_pipeline/mg-process-test1/output/",
        "output_prefix" : "output_test"
        }

    metadata = {
        "chinput_1" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606),
        "chinput_2" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606)
        }

    config = {
        "chicago_setting_file":
            "/Users/pacera/test_pipeline/mg-process-test1/tests/"+
            "data/sGM12878Settings/sGM12878.settingsFile",
        "chicago_desing_dir":
            "/Users/pacera/test_pipeline/mg-process-test1/"+
            "tests/data/hg19TestDesign/",
        #"chicago_print_memory": None,
        "chicago_cutoff": "5",
        "chicago_export_format": "washU_text",
        #"chicago_export_order": None,
        #"chicago_rda": None,
        #"chicago_save_df_only": None,
        "chicago_examples_prox_dist": "1e6",
        #"chicago_examples_full_range": None,
        #"chicago_en_feat_files": None,
        #"chicago_en_feat_folder": None,
        "chicago_en_min_dist": "0",
        "chicago_en_max_dist": "1e6",
        #"chicago_en_full_cis_range": None,
        "chicago_en_sample_no": "100"}
        #"chicago_en_trans": None,
        #"chicago_features_only":None}


    chicago_handle = ChicagoTool(config)
    chicago_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["output_dir"] + "/data/" + output_files["output_prefix"] + ".Rds") is True
    assert os.path.getsize(output_files["output_dir"] + "/data/" + output_files["output_prefix"] + ".Rds") > 0
