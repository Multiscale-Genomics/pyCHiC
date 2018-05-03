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

    def __init__(self, configuration=None):
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

    def makeDesignFiles(self, designDir, outPrefixDesign, parameters):
        """
        make the design files and store it in the specify design folder. It is a
        wrapper of makeDesignFiles.py

        Parameters:
        -----------
        designDir: str,
                   Path to the folder with the output files(recommended the same
                   folder as .map and .baitmap files).
        parameters: dict,
                    list of parameter already selected by
                    get_makeDesignFiles_params().
        outPrefixDesign: str
            Name of the output Design files, recomended
            to be the same name as .rmap and .baitmap files
        Returns:
        -------
        bool

        """
        #if makeDesignFiles.py is added to PATH
        args = ["makeDesignFiles.py", "--outfilePrefix", outPrefixDesign,
            "--designDir", designDir]

        args += parameters
        print(args)

        logger.info("makeDesignFile : "+ " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        if os.path.isfile(outPrefixDesign + ".nbpb") is True:
            pass
        else:
            logger.fatal("makeDesignFiles.py failed to generate design files")
            logger.fatal("makeDesignFiles stderr" + proc_err)
            logger.fatal("makeDesignFiles stdout" + proc_out)
            return False

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
            "makeDesignFiles_rmapfile" : ["--rmapfile", True],
            "makeDesignFiles_baitmapfail" : ["--baitmapFIle", True]
            }

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
            designDir : path to the designDir containin .rmap and .baitmap files
        input_metadata: dict
        output_files: dict
            outFilePrefix : path to the output folder and prefix name of files
                example: "/folder1/folder2/prefixname". Recommended to use the
                path to designDir and the same prefix as .rmap and .baitmap

        Returns:
        --------
        output_files : dict
            List of location for the output files.
        output_metadata : dict
            List of matching metadata dict objects.
        """

        commands_params = self.get_makeDesignFiles_params(self.configuration)

        logger.info("makeDesignFiles command parameters " + " ".join(commands_params))

        results = self.makeDesignFiles(input_files["designDir"],
                                       output_files["outPrefixDesign"],
                                       commands_params)

        results = compss_wait_on(results)

        output_metadata = {
            "output" : Metadata(
                data_type="Designfiles",
                file_type=[".nbpb", ".npb", ".poe"],
                file_path=input_files["designDir"],
                sources=[
                    input_metadata[".rmap"].file_path,
                    input_metadata[".baitmap"].file_path
                    ],
                taxon_id=[
                    input_metadata[".rmap"].taxon_id,
                    input_metadata[".baitmap"].taxon_id
                    ],
                meta_data={
                    "tool" : "makeDesignFiles, make the design files"+
                             "used by chicago as part of the input file"
                }

            )
        }

        return (results, output_metadata)
