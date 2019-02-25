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
from pyCHiC.tool.makeBaitmap import makeBaitmapTool


def test_baitmap():
    """
    Function to test the makeBaitmap tool. Tool that
    generate .baitmap files
    """

    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    path = os.path.join(os.path.dirname(__file__), "data/")

    configuration = {
        "execution": path,
        "chic_RE_name": "HindIII",
        "chic_RE_sequence": "A|AGCTT",
        "bowtie2_fasta_input" : "True"
    }

    input_files = {
        "bowtie_gen_idx": path + "test_baitmap/chr21_hg19.fa.bt2.tar.gz",
        "probes_fa": path + "test_baitmap/h19_promoter.fa",
        "Rtree_file_dat": path + "test_rmap/rtree_file.dat",
        "Rtree_file_idx": path + "test_rmap/rtree_file.idx",
        "genome_fa": path + "test_baitmap/chr21_hg19.fa",
        "chr_handler": path + "test_baitmap/chr_handler.txt"
    }

    output_files = {
        "bait_sam":  path + "test_baitmap/baits.sam",
        "out_bam": path + "test_baitmap/baits.bam",
        "out_baitmap": path + "test_run_chicago/test.baitmap"
    }

    metadata = {
        "bowtie_gen_idx": Metadata(
            "bowtie_gen_idx", "tar", input_files["bowtie_gen_idx"], [input_files["genome_fa"]],
            {
                "assembly": "test",
                "tool": "bowtie_gen_idx"
            }, 9606
            ),
        "genome_fa": Metadata(
            "hg38", "fasta", input_files["genome_fa"], [],
            {
                "assembly": "test",
                "tool": "bowtie_gen_idx"
            }, 9606),

        "probes_fa": Metadata(
            "C-HiC probes", "fasta", input_files["probes_fa"], [],
            {
                "assembly": "test",
                "tool": "bowtie_gen_idx"
            }, 9606),

        "Rtree_file_dat": Metadata(
            "Rtree files", "dat", input_files["Rtree_file_dat"], [],
            {"genome": input_files["genome_fa"],
             "RE": {"HindIII": 'A|AGCTT'}},
            9606
            ),

        "Rtree_file_idx": Metadata(
            "Rtree files", "idx", input_files["Rtree_file_idx"], [],
            {"genome": input_files["genome_fa"],
             "RE": {"HindIII": 'A|AGCTT'}},
            9606
            )
    }

    baitmap_handler = makeBaitmapTool(configuration)
    baitmap_handler.run(input_files, metadata, output_files)

    assert os.path.getsize(output_files["out_bam"]) > 0
    assert os.path.getsize(output_files["out_baitmap"]) > 0
