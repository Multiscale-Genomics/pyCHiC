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
import pytest # pylint: disable=unused-import
import os

from basic_modules.metadata import Metadata
from CHiC.process_bam2chicago_tool import process_bam2chicago

def test_process_bam2chicago():
    """
    Function to test bam2chicago.py
    """
    path = os.path.join(os.path.dirname(__file__),"data/")

    input_files = {
        "RMAP" : path + "test_run_chicago/test.rmap",
        "BAITMAP" : path +  "test_run_chicago/test.baitmap",
        "BAM" : path + "test_bed2bam/outbam_sorted.bam",
    }

    output_files = {
        "sample_name" :  path + "test_bam2chicago_tool/sampleout",
        "chrRMAP" : path + "test_bam2chicago_tool/chrtest.rmap",
        "chrBAITMAP" : path + "test_bam2chicago_tool/chrtest.baitmap"
    }

    metadata = {
        "RMAP" : Metadata(
            "data_chicago_input", ".rmap",
            path+"/h19_chr20and21_chr.rmap", None, {}, 9606),
        "BAITMAP" : Metadata(
            "data_chicago_input", ".baitmap",
            path+"/h19_chr20and21.baitmap_4col_chr.txt", None, {}, 9606),
        "BAM" : Metadata(
            "txt", "bamfile", path + "/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606)
    }

    bam2chicago_handle = process_bam2chicago()
    bam2chicago_handle.run(input_files, metadata, output_files)

    out_path = output_files["sample_name"] + "/sampleout.chinput"

    assert os.path.isfile(out_path) is True
    assert os.path.getsize(out_path) > 0
