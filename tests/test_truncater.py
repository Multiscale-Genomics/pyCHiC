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

from tool.truncater import Truncater

def test_truncater():
    """
    Test for the truncater function
    """

    path = os.path.join(os.getcwd(),"data/test_truncater/")

    input_files = {
        "fastq1" : path + "SRR3535023_1.fastq",
        "fastq2" : path + "SRR3535023_2.fastq"
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
            )

        }

    output_files = {
        "outdir" : path
    }

    configuration = {
        "quiet_progress": True,
        "RE": "A^AGCT",
        "threads" : "2"
        }

    truncater_hdl = Truncater(configuration)
    truncater_hdl.run(input_files, input_metadata, output_files)

    assert os.path.isfile(output_files["outdir"]+"SRR3535023_1.trunc.fastq") is True
    assert os.path.isfile(output_files["outdir"]+"SRR3535023_2.trunc.fastq") is True

    assert os.path.getsize(output_files["outdir"]+"SRR3535023_1.trunc.fastq") > 0
    assert os.path.getsize(output_files["outdir"]+"SRR3535023_2.trunc.fastq") > 0
