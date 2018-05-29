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
from basic_modules.metadata import Metadata
import os

from process_chicago_CHiC import process_chicago_CHiC

def test_process_CHiC():
    """
    Test for the process_chicago_CHiC
    """
    path = os.path.join(os.getcwd(),"data/test_process_CHiC/")

    input_files = {
        "fastq1" : path + "SRR3535023_1.fastq",
        "fastq2" : path + "SRR3535023_2.fastq",
        "genome_fa" : path + "cc"
    }

    input_metadata = {
            "fastq1": Metadata(
                data_type="text",
                file_type="fastq",
                file_path=input_files["fastq1"],
                sources="",
                taxon_id=9606,
                meta_data=""
            ),
            "fastq2": Metadata(
                data_type="text",
                file_type="fastq",
                file_path=input_files["fastq2"],
                sources="",
                taxon_id=9606,
                meta_data=""
            ),
            "genome_fa" : Metadata(
                data_type="text",
                file_type="fasta",
                file_path=input_files["genome_fa"],
                sources="",
                taxon_id=9606,
                meta_data="GRCh38",
            )
        }

    output_files = {
        "out_dir" : path + "output/",
        "out_prefix_makeRmap" : "restriction_enzyme_test_HindIII_hg38.txt",
        "Rtree_files" : path + "output/rtree_file"
    }

    configuration = {
        "RE_truncater": "A^AGCT",
        "RE": {"HindIII" : 'A|AGCTT'}
        }

    CHiC_hdl = process_chicago_CHiC(configuration)
    CHiC_hdl.run(input_files, input_metadata, output_files)

    #assert truncater.py

    assert os.path.isfile(output_files["out_dir"]+"SRR3535023_1.trunc.fastq") is True
    assert os.path.isfile(output_files["out_dir"]+"SRR3535023_2.trunc.fastq") is True

    assert os.path.getsize(output_files["out_dir"]+"SRR3535023_1.trunc.fastq") > 0
    assert os.path.getsize(output_files["out_dir"]+"SRR3535023_2.trunc.fastq") > 0

    #assert makeRmap_Tool.py

    out = "".join(
        [
            f for f in os.listdir(output_files["out_dir"])
            if f.startswith("Digest_") and f.endswith(".map")
        ]
    )

    assert os.path.getsize(output_files["out_dir"] + out) > 0
    assert os.path.getsize(output_files["Rtree_files"] + ".dat")
    assert os.path.getsize(output_files["Rtree_files"] + ".idx")

    #assert makeBaitmap.py

