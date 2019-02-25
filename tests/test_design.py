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
import os.path
import pytest # pylint: disable=unused-import


from basic_modules.metadata import Metadata

from pyCHiC.tool.makeDesignFiles import makeDesignFilesTool

def test_design():
    """
    Function to test the right generation of CHicago Design files
    """
    path = os.path.join(os.path.dirname(__file__), "data/")


    config_file = {
        "makeDesignFiles_minFragLen" : "150",
        "makeDesignFiles_maxFragLen" : "40000",
        "makeDesignFiles_maxLBrownEst" : "1500000",
        "makeDesignFiles_binsize" : "20000",
        "makeDesignFiles_removeb2b" : True,
        "makeDesignFiles_removeAdjacent" : True,
        "makeDesignFiles_outfilePrefix" : path + "test_run_chicago/test",
        #"makeDesignFiles_designDir" : path + "test_run_chicago",
        "makeDesignFiles_rmap" : path + "test_run_chicago/test.rmap",
        "makeDesignFiles_baitmap" :  path + "test_run_chicago/test.baitmap"
        }

    input_files = {
        "RMAP" : path + "test_run_chicago/test.rmap",
        "BAITMAP": path + "test_run_chicago/test.baitmap"
    }

    metadata = {
        "RMAP" : Metadata(
            "data_chicago_input", ".rmap",
            path + "test_run_chicago", None, {}, 9606),
        "BAITMAP" : Metadata(
            "data_chicago_input", ".baitmap",
            path + "test_run_chicago", None, {}, 9606)
    }

    output_files = {
        "nbpb" : path + "test_run_chicago/test.nbpb",
        "npb" : path + "test_run_chicago/test.npb",
        "poe" : path + "test_run_chicago/test.poe"
    }

    design_handle = makeDesignFilesTool(config_file)
    design_handle.run(input_files, metadata, output_files)

    assert os.path.isfile(path + "test_run_chicago/test" + ".nbpb") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".nbpb") > 0

    assert os.path.isfile(path + "test_run_chicago/test" + ".npb") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".npb") > 0

    assert os.path.isfile(path + "test_run_chicago/test" + ".poe") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".poe") > 0
