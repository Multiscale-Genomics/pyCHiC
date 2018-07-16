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

import os.path
import pytest  # pylint: disable=unused-import

from basic_modules.metadata import Metadata
from process_run_chicago import process_run_chicago

def test_process_chicago():
    """
    Test the chicago Pipeline

    Running the chicago pipeline with the test data from the command line
    """


    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "chinput": path + "test_run_chicago/data_chicago/GM_rep1.chinput",
        "setting_file" : path + "test_run_chicago/data_chicago/sGM12878.settingsFile"
        }

    output_files = {
        "output": path + "test_run_chicago/data_chicago/out_run_chicago.tar",
        }

    metadata = {
        "chinput" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606)
        }

    config = {
        "chicago_design_dir": path + "/test_run_chicago/data_chicago",
        "chicago_print_memory": "None",
        "chicago_out_prefix" : "output_test",
        "chicago_cutoff": "5",
        "chicago_export_format": "washU_text",
        "chicago_export_order": "None",
        "chicago_rda": "None",
        "chicago_save_df_only": "None",
        "chicago_examples_prox_dist": "1e6",
        "chicago_examples_full_range": "None",
        "chicago_en_feat_files": "None",
        "chicago_en_min_dist": "0",
        "chicago_en_max_dist": "1e6",
        "chicago_en_full_cis_range": "None",
        "chicago_en_sample_no": "100",
        "chicago_en_trans": "None",
        "chicago_features_only": "None"}


    chicago_handle = process_run_chicago(config)
    chicago_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["output"]) is True

    assert os.path.getsize(output_files["output"]) > 0

