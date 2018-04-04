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
from tool.makeRmapFileTool import makeRmapFile


def test_makeRmapFileTool():
    """
    Function to test generation of .rmap input files
    from CHiCAGO
    """

    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "genome" : path + "test_makeRmap/mm10.fa"
    }

    output_files = {
    "output_dir" : path + "test_makeRmap/",
    "output_prefix" : "test_digest"
    }

    metadata = {
        "genome" : Metadata(
        "hg38", "fasta", "../genome_mm10/mm10.fa", None,
        "HindIII", 9606),
    }

    config = {
        "digester_re1" : "A^AGCTT,HindIII"
    }

    makeRmap_handle = makeRmapFile(config)
    makeRmap_handle.run(input_files, metadata, output_files)

    out =  "".join([f for f in os.listdir(output_files["output_dir"])
            if f.startswith("Digest_") and f.endswith(".map")])

    assert os.path.getsize(output_files["output_dir"] + out) > 0
