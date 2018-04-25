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
import os.path

from basic_modules.metadata import Metadata
from tool.bam2chicago_Tool import bam2chicagoTool


def test_bam2chicagoTool():
    """
    Function to test if bam2chicago convert bam files to .chinput files
    """

    path = os.path.join(os.path.dirname(__file__), "data/")

    input_metadata = {
        "bam_1": Metadata(
            "input_chicago", "chinput", [], None, {"assembly" : "test"}, 9606),
        "bsgenome" : Metadata(
            "data_bam2chicago", "bsgenome", [], None, {"assembly" : "test"}, 9606)
        }

    input_files = {
        "bamFile" : path + "test_bam2chicago/SRR1658573_merge.bamx",
        "rmapFile" : path + "test_Design/h19_chr20and21.rmap",
        "baitmapFile" : path + "test_Design/h19_chr20and21.baitmap"
    }

    output_files = {
        "sample_name" : "/Users/pacera/developing/C-HiC/tests/data/test_bam2chicago/test_bam2chicago.chinput"
    }


    test1 = bam2chicagoTool()
    test1.run(input_files, input_metadata, output_files)

    #write the assertions
