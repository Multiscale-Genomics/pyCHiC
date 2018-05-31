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
from tool.makeRmap_tool import makeRmapFile


def test_makeRmapFileTool():
    """
    Function to test generation of .rmap input files
    from CHiCAGO
    """
    path = os.path.join(os.path.dirname(__file__), "data/")


    configuration = {"RE" : {"HindIII" : 'A|AGCTT'}
                    }

    input_files = {
        "genome_fa" : path+ "test_makeBaitmap/toy_GRCh38.fa",
        }


    metadata = {"genome_fa" : Metadata(
        "txt", "fasta", path+ "test_makeBaitmap/toy_GRCh38.fa",
        None, 9606, ""),
        }


    output_files = {
        "out_dir_makeRmap" : path + "test_Design/",
        "out_prefix_makeRmap" : "test",
        "Rtree_files" : path + "test_makeRmap/rtree_file"
        }


    makeRmap_handle = makeRmapFile(configuration)
    makeRmap_handle.run(input_files, metadata, output_files)

    out = "".join(
        [
            f for f in os.listdir(output_files["out_dir_makeRmap"])
            if f.startswith("Digest_") and f.endswith(".map")
        ]
    )

    assert os.path.getsize(output_files["out_dir_makeRmap"] + out) > 0
    assert os.path.getsize(output_files["Rtree_files"] + ".dat")
    assert os.path.getsize(output_files["Rtree_files"] + ".idx")
