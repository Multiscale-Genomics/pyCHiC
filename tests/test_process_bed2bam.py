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

from process_bed2bam import process_bed2bam
from basic_modules.metadata import Metadata

def test_process_bed2bam():
    """
    Test for bed2chicagobamWra.py
    """

    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "bed" : path + "test_fastq2bed/03_filtered_reads/valid_r1-r2_intersection_b51cdf1282.tsv",
        "ncpus" : "2"
    }

    input_metadata = {
            "bed": Metadata(
                data_type="text",
                file_type="tsv",
                file_path=input_files["bed"],
                sources="",
                taxon_id=9606,
                meta_data=""
            )
        }

    output_files = {
        "bam_out" : path + "test_bed2bam/outbam"
    }

    bed2bam_hdl = process_bed2bam()
    bed2bam_hdl.run(input_files, input_metadata, output_files)

    assert os.path.isfile(output_files["bam_out"]+"_sorted.bam") is True
    assert os.path.isfile(output_files["bam_out"]+"_sorted.bam") is True
