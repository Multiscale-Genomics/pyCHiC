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
import os
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata
from tool.makeBaitmap import makeBaitmapTool

def test_baitmap():
    """
    Function to test the makeBaitmap tool. Tool that
    generate .baitmap files
    """
    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "genome_idx" : path + "test_baitmap/chr21_hg19.fa",
        "probes_fa" : path + "test_baitmap/baits.fa",
        "Rtree_files" : path + "test_rmap/rtree_file"
    }

    output_files = {
        "out_sam" :  path + "test_baitmap/baits.sam",
        "out_baitmap" : path + "test_run_chicago/test.baitmap"
    }

    metadata = {
        "genome_digest" : Metadata(
            "hg38", "fasta", path + "test_rmap/chr21_hg19.fa",
            None, "HindIII", 9606),

        "probes" : Metadata(
            "C-HiC probes", "fasta", path + "test_baitmap/baits.fa",
            None, None, 9606),

        "Rtree_files" : Metadata(
            "Rtree files", [".dat", ".idx"], path + "test_rmap/rtree_file",
            {"genome" : path + "test_rmap/chr21_hg19.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606
            )
    }

    makeBaitmap_handler = makeBaitmapTool()
    makeBaitmap_handler.run(input_files, metadata, output_files)

    assert os.path.getsize(output_files["out_sam"]) > 0
    assert os.path.getsize(output_files["out_baitmap"]) > 0
