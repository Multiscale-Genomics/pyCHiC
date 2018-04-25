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

from process_baitmap import generate_CHiCAGO_baitmap
from basic_modules.metadata import Metadata

def test_process_rmapBaitmap():
    """
    Test for process_rmapBaitmap pipeline.
    This pipeline generate .baitmap files,
    input files for CHiCAGO pipeline
    """

    path = os.path.join(os.path.dirname(__file__) + "/data")

    configuration = {"RE" : {"HindIII" : 'A|AGCTT'},
                    }

    input_files = {
        "genome" :  path + "/test_makeBaitmap/toy_GRCh38.fa",
        "probes_fa": path + "/test_makeBaitmap/baits.fa",
        "Rtree_files" : path + "/test_process_rmap/rtree_file",
        }

    input_metadata = {
        "probes" : Metadata(
            "C-HiC probes", "fasta", path + "/test_makeBaitmap/baits.fa",
            None, None, 9606),

        "Rtree_files" : Metadata(
            "Rtree files", [".dat", ".idx"], path + "/test_makeRmap/rtree_file",
            {"genome" : path + "/test_makeRmap/toy_GRCh38.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606),

        "genome_digest" : Metadata(
            "hg38", "fasta", path + "/test_makeRmap/toy_GRCh38.fa", None, "HindIII", 9606),
        }

    output_files = {
        "out_sam" :  path + "/test_process_baitmap/baits.sam",
        "out_baitmap" : path + "/test_process_baitmap/test.baitmap"
    }

    generate_CHiCAGO_baitmap_hand = generate_CHiCAGO_baitmap(configuration)
    generate_CHiCAGO_baitmap_hand.run(input_files, input_metadata, output_files)

    assert os.path.getsize(output_files["out_baitmap"]) > 0
