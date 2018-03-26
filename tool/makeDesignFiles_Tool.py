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
import subprocess
import sys
from utils import logger
import shutil

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

#################################################################

class makeDesignFilesTool(Tool):
    """
    Tool for makeing the design files as part of the input for Chicago
    capture Hi-C
    """

    def __init__(self, configuration={}):
        """
        Initialise the tool with the configuration file

        Parameters:
        -----------
        configuration: dict
            Dictionary containing parameters defining how the tool
            should work
        """
        print("Initialising makeDesingFiles")

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def makeDesignFiles(self, rmapFile, baitMapFile, outFilePrefix,
                        designDir, parameters):
        """
        make the design files and store it in the specify design folder.

        Parameters:
        -----------
        rmapFIle :  A tab-separated file of the format
                    <chr> <start> <end> <numeric ID>,
                    describing the restriction digest (or "virtual digest"
                    if pooled fragments are used). These numeric IDs are referred to as
                    "otherEndID" in Chicago. All fragments mapping outside of the digest
                    coordinates will be disregarded by both these scripts and Chicago.
        baitMapFile: Tab-separated file of the format
                     <chr> <start> <end> <numeric ID> <annotation>,
                     listing the coordinates of the baited/captured
                     restriction fragments (should be a subset of the fragments
                     listed in rmapfile), their numeric IDs (should match those listed
                     in rmapfile for the corresponding fragments) and their annotations
                     (such as, for example, the names of baited promoters). The numeric
                     IDs are referred to as "baitID" in Chicago.
        outFilePrefix: Prefix name of the output files
        designDir: Path to the folder with the output files(recommended the same
                    folder as .map and .baitmap files).
        parameters: list of parameter already selected by
                    get_makeDesignFiles_params().
        Returns:
        -------
        bool
            writes the output files in the defined location

        """
        #if makeDesignFiles.py is added to PATH
        args = ["makeDesignFiles.py", "--rmapfile", rmapFile, "--baitmapfile", baitMapFile,
                "--outfilePrefix", outFilePrefix, "--designDir", designDir]
        args += parameters

        logger.info("makeDesignFile : "+ " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        if os.path.isfile(outFilePrefix + ".nbpb") is True:
            pass
        else:
            logger.fatal("makeDesignFiles.py failed to generate design files")
            logger.fatal("makeDesignFiles stderr" + proc_err)
            logger.fatal("makeDesignFiles stdout" + proc_out)
            return False

        #if design
        return True

    @staticmethod
    def get_makeDesignFiles_params(params):
        """
        This function handle chicago parameters,
        selecting the given ones and passing to the
        command line.
        """

        command_params = []

        command_parameters = {
            "makeDesignFiles_minFragLen" : ["--minFragLen", True],
            "makeDesignFiles_maxFragLen" : ["--maxFragLen", True],
            "makeDesignFiles_maxLBrownEst" :["--maxLBrownEst", True],
            "makeDesignFiles_binSize" : ["--binSize", True],
            "makeDesignFiles_removeb2b" : ["--removeb2b", False],
            "makeDesignFiles_removeAdjacent" : ["--removeAdjacent", False],
            }
        """
        "makeDesignFiles_rmapFile" : ["-r", True],
        "makeDesignFiles_baitMapFile" : ["-f", True],
        "makeDesignFiles_outFilePrefix" : ["-o", True],
        "makeDesignFiles_designDir" : ["-d", True]
        """


        for parameter in params:
            if parameter in command_parameters:
                if command_parameters[parameter][1]:
                    command_params += [command_parameters[parameter][0], params[parameter]]
                else:
                    command_params += [command_parameters[parameter][0]]

        print(command_params)
        return command_params


    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run makeDesignFiles.

        Parameters:
        ----------

        input_files: dict
            rmapfile : location
            baitMapFile : location
        input_metadata: dict
        output_files: dict
            outFilePrefix : prefix output name
            designDir : location of output folder

        Returns:
        --------
        output_files : dict
            List of location for the output files.
        output_metadata : dict
            List of matching metadata dict objects.
        """

        if not os.path.exists(output_files["designDir"]):
            logger.error(output_files["designDir"] + "does not" +
                + "exists")
            logger.info("Introduce a valid directory with .map and .mapbait")

        commands_params = self.get_makeDesignFiles_params(self.configuration)

        logger.info("makeDesignFiles command parameters " + " ".join(commands_params))

        results = self.makeDesignFiles(input_files["rmapFile"],
                                 input_files["baitMapFile"],
                                  output_files["outFilePrefix"],
                                  output_files["designDir"],
                                  commands_params)

        output_metadata ={
            "output" : Metadata(
                data_type = "Designfiles",
                file_type = [".nbpb", ".npb", ".poe"],
                file_path = output_files["designDir"],
                sources = [
                    input_metadata[".rmap"].file_path,
                    input_metadata[".baitmap"].file_path
                    ],
                taxon_id = [
                    input_metadata[".rmap"].taxon_id,
                    input_metadata[".baitmap"].taxon_id
                    ],
                meta_data = {
                    "tool" : "makeDesignFiles, make the design files"+
                        "used by chicago as part of the input file"
                }

            )
        }

        return (results, output_metadata)




config_file = {
    "makeDesignFiles_minFragLen" : "150",
    "makeDesignFiles_maxFragLen" : "40000",
    "makeDesignFiles_maxLBrownEst" : "1.5e6",
    "makeDesignFiles_binSize" : "20000",
    "makeDesignFiles_removeb2b" : True,
    "makeDesignFiles_removeAdjacent" : True,
    }

input_files = {
    "rmapFile" : "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/PCHiCdata/inst/extdata/hg19TestDesign/h19_chr20and21.rmap",
    "baitMapFile" : "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/PCHiCdata/inst/extdata/hg19TestDesign/h19_chr20and21.baitmap"
}

input_metadata = {
    ".rmap" : Metadata(
        "data_chicago_input", ".rmap",
         "/Users/pacera/test_makedir/output/h19_chr20and21.rmap",
         None, {}, 9606),
    ".baitmap" : Metadata(
        "data_chicago_input", ".rmap",
         "/Users/pacera/test_makedir/output/h19_chr20and21.baitmap",
         None, {}, 9606)
}

output_files = {
    "outFilePrefix" : "/Users/pacera/designTest/testecillo",
    "designDir" : "/Users/pacera/designTest"
}

test = makeDesignFilesTool(config_file)

print(test.run(input_files, input_metadata, output_files))



































