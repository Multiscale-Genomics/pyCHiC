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
from tool.bam2chicago import bam2chicago

def test_bam2chicago():
	"""
	Function to test bam2chicago.py
	"""
	path = "data/test_bam2chicago"

	input_files = {
		"RMAP" : path + "/h19_chr20and21_chr.rmap",
		"BAITMAP" : path +  "/h19_chr20and21.baitmap_4col_chr.txt",
		"BAM" : path + "/SRR3535023_1_2.hicup.bam"
	}

	output_files = {
		"out_dir" : path,
		"output_file" :  "/sample"
	}

	input_metadata = {}

	bam2chicago_handle = bam2chicago()
	bam2chicago_handle.run(input_files, input_metadata, output_files)