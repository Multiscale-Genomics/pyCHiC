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
import os
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata
from CHiC.tool.hicup_tool import hicup

def test_hicup():
    """
    This function is goin to test that hicup works and produce the right output
    """
    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "genome_fa" : path + "test_baitmap/chr21_hg19.fa",

        "fastq1" : path + "test_truncater/SRR3535023_1.fastq",
        "fastq2" : path + "test_truncater/SRR3535023_2.fastq",
        "bowtie_gen_idx" : path + "test_baitmap/chr21_hg19.fa.bt2.tar.gz"
    }

    configuration = {
        "hicup_renzyme" : "A^AGCTT,HindIII",
        "genome_name" : "test_hg19",
        "hicup_bowtie2_loc": "/usr/bin/bowtie2",
        "hicup_longest": "800",
        "hicup_shortest": "150",
        "hicup_outdir": path + "test_hicup/output",
        "hicup_zip": "True",
    }

    output_files = {
        "hicup_outdir_tar" : path + "test_hicup/output.tar"
    }

    metadata = {
        "genome_fa" : Metadata(
            "TXT", "FASTA",
            input_files["genome_fa"], None, {}, 9606),

        "fastq1" : Metadata(
            "TXT", "FASTQ",
            input_files["fastq1"], None, {}, 9606),


        "fastq2" : Metadata(
            "TXT", "FASTQ",
            input_files["fastq2"], None, {}, 9606)
    }

    hicup_handl = hicup(configuration)
    hicup_handl.run(input_files, metadata, output_files)

    assert os.path.isfile(output_files["hicup_outdir_tar"]) is True
    assert os.path.getsize(output_files["hicup_outdir_tar"]) > 0
