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
import pytest

from process_rmap import generate_CHiCAGO_rmap
from basic_modules.metadata import Metadata

def test_process_rmap():
    """
    Test for process_rmapBaitmap pipeline.
    This pipeline generate .rmap file,
    input files for CHiCAGO pipeline
    """

    path = os.path.join(os.path.dirname(__file__) + "/data")

    configuration = {"RE" : {"HindIII" : 'A|AGCTT'},
                    }

    input_files = {
        "genome":  path + "/test_makeRmap/toy_GRCh38.fa",
        }

    input_metadata = {

        "Rtree_files" : Metadata(
            "Rtree files", [".dat", ".idx"], path + "/test_makeRmap/rtree_file",
            {"genome" : path + "/test_makeRmap/toy_GRCh38.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606),

        "genome_digest" : Metadata(
            "hg19", "fasta", path + "/test_makeRmap/toy_GRCh38.fa", None, "HindIII", 9606),
        }

    output_files = {
        "out_dir_makeRmap" : path + "/test_process_rmap/",
        "out_prefix_makeRmap" : "restriction_enzyme_test2",
        "Rtree_files" : path + "/test_process_rmap/rtree_file",
    }

    generate_CHiCAGO_rmap_hand = generate_CHiCAGO_rmap(configuration)
    generate_CHiCAGO_rmap_hand.run(input_files, input_metadata, output_files)

    assert os.path.getsize(output_files["out_dir_makeRmap"] +
                           output_files["out_prefix_makeRmap"] + ".rmap") > 0
