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
import gzip
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool.gem_indexer import gemIndexerTool

@pytest.mark.hic
@pytest.mark.genome
def test_gem_indexer():
    """
    Test case to ensure that the GEM indexer works.
    """

    genome_gem_idx = "/home/pablo/MuG/C-HiC/tests/data/test_gem_indexer/tb.Human.GCA_000001405.21.fasta.gem"

    input_files = {
        "genome": "/home/pablo/MuG/data/chr21.fa"
    }

    output_files = {
        "index": "/home/pablo/MuG/C-HiC/tests/data/test_gem_indexer/tb.Human.GCA_000001405.21.fasta.gem.gz"
    }

    metadata = {
        "genome": Metadata(
            "Assembly", "fasta", None, 9606,
            {'assembly' : 'test'},
            )
    }


    gem_it = gemIndexerTool()
    gem_it.run(input_files, metadata, output_files)

    assert os.path.isfile(genome_gem_idx) is True
    assert os.path.getsize(genome_gem_idx) > 0
