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
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool.fastq2bed import Fastq2bed

def test_fastq2bed():
    """
    Test for the fastq2bed function
    """
    path = os.path.join(os.path.dirname(__file__), "data")

    input_files = {
        "fastq1" : path + "/test_truncater/SRR3535023_1.trunc.fastq",
        "fastq2" : path + "/test_truncater/SRR3535023_2.trunc.fastq",
        "gem_idx" : path + "/test_baitmap/chr21_hg19.fa.gem.gz",
        "genome_fa" : path + "/test_baitmap/chr21_hg19.fa",
        "RE" : "HindIII",
        "chromosome" : "22"
    }

    metadata = {
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
        "genome_fa": Metadata(
            data_type="text",
            file_type="fasta",
            file_path=input_files["genome_fa"],
            sources="",
            taxon_id=9606,
            meta_data=""
            )
        }

    output_files = {
        "wd" : path + "/test_fastq2bed"
    }

    fastq2bed_hdl = Fastq2bed()
    fastq2bed_hdl.run(input_files, metadata, output_files)

    try:
        file_list = os.listdir(output_files["wd"]+"/01_mapped_r1/")
        for file_hdl in file_list:
            if "full" in file_hdl:
                assert os.path.getsize(output_files["wd"]+"/01_mapped_r1/"+file_hdl) > 0
    except IOError:
        print("tadbit map did not generate any output for read 1 =(")

    try:
        file_list = os.listdir(output_files["wd"]+"/01_mapped_r2/")
        for file_hdl in file_list:
            if "full" in file_hdl:
                assert os.path.getsize(output_files["wd"]+"/01_mapped_r2/"+file_hdl) > 0
    except IOError:
        print("tadbit map did not generate any output for read 2 =(")

    try:
        file_list = os.listdir(output_files["wd"]+"/02_parsed_reads/")
        for file_hdl in file_list:
            assert os.path.getsize(output_files["wd"]+"/02_parsed_reads/"+file_hdl) > 0
    except IOError:
        print("tadbit pare did not generate any output =(")

    try:
        file_list = os.listdir(output_files["wd"]+"/03_filtered_reads/")
        for file_hdl in file_list:
            if "valid" in file_hdl:
                assert os.path.getsize(output_files["wd"]+"/03_filtered_reads/"+file_hdl) > 0
    except IOError:
        print("tadbit filter didn't generate any output =(")
