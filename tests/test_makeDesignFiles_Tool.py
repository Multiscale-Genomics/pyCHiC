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

from tool.makeDesignFiles_Tool import makeDesignFilesTool

def test_makeDesignFilesTool():
    """
    Function to test the right generation of CHicago Design files
    """
    path = os.path.join(os.path.dirname(__file__), "data/test_Design")


    config_file = {
        "makeDesignFiles_minFragLen" : "150",
        "makeDesignFiles_maxFragLen" : "40000",
        "makeDesignFiles_maxLBrownEst" : "1.5e6",
        "makeDesignFiles_binSize" : "20000",
        "makeDesignFiles_removeb2b" : True,
        "makeDesignFiles_removeAdjacent" : True,
        #"rmapFile" : "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/"+
        #             "PCHiCdata/inst/extdata/hg19TestDesign/h19_chr20and21.rmap",
        #"baitMapFile" : "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/"+
        #                "PCHiCdata/inst/extdata/hg19TestDesign/h19_chr20and21.baitmap"
        }

    input_files = {
        "designDir" : path
    }

    input_metadata = {
        ".rmap" : Metadata(
            "data_chicago_input", ".rmap",
            path, None, {}, 9606),
        ".baitmap" : Metadata(
            "data_chicago_input", ".baitmap",
            path, None, {}, 9606)
    }

    output_files = {
        "outPrefixDesign" : path + "/h19_chr20and21_test",
    }

    makeDesignFiles_handle = makeDesignFilesTool(config_file)
    makeDesignFiles_handle.run(input_files, input_metadata, output_files)

    assert os.path.isfile(path + "/h19_chr20and21_test" + ".nbpb") is True
    assert os.path.getsize(path + "/h19_chr20and21_test" + ".nbpb") > 0

    assert os.path.isfile(path + "/h19_chr20and21_test" + ".npb") is True
    assert os.path.getsize(path + "/h19_chr20and21_test" + ".npb") > 0

    assert os.path.isfile(path + "/h19_chr20and21_test" + ".poe") is True
    assert os.path.getsize(path + "/h19_chr20and21_test" + ".poe") > 0
