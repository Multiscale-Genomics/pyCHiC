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
import argparse
import sys

import pytest  # pylint: disable=unused-import


def all_toolchain(verbose=False):
    """
    Runs the tests for all of the tools

    This set is only required for determining code coverage.
    run from home of the repo
    """
    directories = ["./data/test_bam2chicago_tool",
                   "./data/test_bed2bam",
                   "./data/test_design",
                   "./data/test_fastq2bed",
                   "./data/test_rmap"]

    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

    params = []

    if verbose is True:
        params.append('-')

   #params.append('test_gem_indexer.py')
    params.append('test_truncater.py')
    params.append('test_rmap_tool.py')
    params.append('test_baitmap.py')
    params.append('test_design.py')
    params.append('test_fastq2bed.py')
    params.append('test_bed2bamchicago.py')
    params.append('test_bam2chicago_tool.py')
    params.append('test_run_chicago.py')

    return pytest.main(params)


if __name__ == "__main__":

    sys._run_from_cmdl = True # pylint: disable=protected-access

    PARSER = argparse.ArgumentParser(description="Test for running all tools in a chain")
    PARSER.add_argument("--verbose", action="store_const", const=True, default=False)

    ARGS = PARSER.parse_args()

    VERBOSE = ARGS.verbose

    all_toolchain()
