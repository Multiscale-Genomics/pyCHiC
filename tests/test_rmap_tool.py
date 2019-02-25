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
import pytest  # pylint: disable=unused-import

from basic_modules.metadata import Metadata
from pyCHiC.tool.rmap_tool import makeRmapFile


def test_rmap_tool():
    """
    Function to test generation of .rmap input files
    from CHiCAGO
    """
    path = os.path.join(os.path.dirname(__file__), "data/")

    configuration = {"chic_RE_name" : "HindIII",
                     "chic_RE_sequence" : "A|AGCTT"
                    }

    input_files = {
        "genome_fa": path + "test_baitmap/chr21_hg19.fa",
        }

    metadata = {
        "genome_fa": Metadata(
            "txt", "fasta", path + "test_baitmap/chr21_hg19.fa",
            None, 9606, ""),
    }

    output_files = {
        "RMAP": path + "test_run_chicago/test.rmap",
        "Rtree_file_dat": path + "test_rmap/rtree_file.dat",
        "Rtree_file_idx": path + "test_rmap/rtree_file.idx",
        "chr_handler": path + "test_baitmap/chr_handler.txt"
        }

    rmap_handle = makeRmapFile(configuration)
    rmap_handle.run(input_files, metadata, output_files)

    assert os.path.getsize(output_files["Rtree_file_dat"]) > 0
    assert os.path.getsize(output_files["Rtree_file_idx"]) > 0
