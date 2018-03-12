"""
Lisence

test the tool practical_tool_pablol.py
"""

from __future__ import print_function

import os.path

from practical_tool_pablo import bwaIndexerTool

def test_tool():
    """
    function to test if the bwa indexer works
    """

    resource_path = os.path.dirname(__file__)

    text_file = resource_path + "/output.tar.gz"

    input_files = {
    "genome" : "/Users/pacera/reference_genome/chr21.fa"
    }

    output_files= {
    "index" : text_file}

    metadata = {}

    print(input_files, output_files)


    tt_handle = bwaIndexerTool()
    tt_files, tt_meta = tt_handle.run(input_files, metadata, output_files)

    assert output_files["index"] == tt_files["index"]
    print(os.path.isfile(text_file))
    assert os.path.isfile(text_file) is True
    assert os.path.getsize(text_file) > 0