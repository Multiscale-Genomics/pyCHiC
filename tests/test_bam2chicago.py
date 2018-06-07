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
import pytest
import os

from basic_modules.metadata import Metadata
from tool.bam2chicago import bam2chicago

def test_bam2chicago():
    """
    Function to test bam2chicago.py
    """
    path = os.path.join(os.path.dirname(__file__),"data/")

    input_files = {
        "RMAP" : path + "test_bam2chicago_Tool/chrtest.rmap",
        "BAITMAP" : path +  "test_bam2chicago_Tool/chrtest.baitmap",
        "BAM" : path + "test_bed2bam/outbam_sorted.bam"
    }

    output_files = {
        "out_dir" : path,
        "output_file" :  "../scripts/out_py"
    }

    input_metadata = {
        ".rmap" : Metadata(
            "data_chicago_input", ".rmap",
            path+"/h19_chr20and21_chr.rmap", None, {}, 9606),
        ".baitmap" : Metadata(
            "data_chicago_input", ".baitmap",
            path+"/h19_chr20and21.baitmap_4col_chr.txt", None, {}, 9606),
        "bam" : Metadata(
            "txt", "bamfile", path + "/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606)
    }

    bam2chicago_handle = bam2chicago()
    bam2chicago_handle.run(input_files, input_metadata, output_files)
