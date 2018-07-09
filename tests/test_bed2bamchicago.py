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

from basic_modules.metadata import Metadata

from CHiC.tool.bed2bam import bed2bam

def test_bed2bam():
    """
    Test for bed2chicagobamWra.py
    """

    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "bed" : path + "test_fastq2bed/03_filtered_reads/valid_r1-r2_intersection_b51cdf1282.tsv"
    }

    metadata = {
        "bed": Metadata(
            data_type="TXT",
            file_type="tsv",
            file_path=input_files["bed"],
            sources="",
            taxon_id=9606,
            meta_data={
                "visible": True,
                "validated": 1
            }
        )
    }

    output_files = {
        "bam_out" : path + "test_bed2bam/outbam.bam",
        "bam_out_sorted": path + "test_bed2bam/outbam_sorted.bam"
    }

    config_file = {"ncpus" : "2"}

    bed2bam_hdl = bed2bam(config_file)
    bed2bam_hdl.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["bam_out_sorted"]) is True
    assert os.path.getsize(output_files["bam_out_sorted"]) > 0
