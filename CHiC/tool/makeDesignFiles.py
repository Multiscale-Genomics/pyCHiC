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
        Tool.__init__(self)
        print("Initialising makeDesingFiles")

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, RMAP=FILE_IN, BAITMAP=FILE_IN, nbpb=FILE_OUT,
          npb=FILE_OUT, poe=FILE_OUT, parameters=IN)
    def makeDesignFiles(self, RMAP, BAITMAP, nbpb, npb, poe, parameters):
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
        Returns:
        -------
        bool
        outFilePrefix: str
            writes the output files in the defined location

        """
        logger.info("oooooloooocooooooo")
        script = os.path.join(os.path.dirname(__file__), "scripts/makeDesignFiles.py")

        args = ["python", script]

        args += parameters

        logger.info("makeDesignFile : "+ " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        return True

    @staticmethod
    def get_design_params(params):
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
            "makeDesignFiles_binsize" : ["--binsize", True],
            "makeDesignFiles_removeb2b" : ["--removeb2b", False],
            "makeDesignFiles_removeAdjacent" : ["--removeAdjacent", False],
            "makeDesignFiles_rmapfile" : ["--rmapfile", True],
            "makeDesignFiles_baitmapfail" : ["--baitmapFIle", True],
            "makeDesignFiles_outfilePrefix" : ["--outfilePrefix", True],
            "makeDesignFiles_designDir" : ["--designDir", True],
            "makeDesignFiles_rmap" : ["--rmapfile", True],
            "makeDesignFiles_baitmap": ["--baitmapfile", True]
            }

        for parameter in params:
            if parameter in command_parameters:
                if command_parameters[parameter][1]:
                    command_params += [command_parameters[parameter][0], params[parameter]]
                else:
                    command_params += [command_parameters[parameter][0]]

        #print(command_params)
        return command_params

    def run(self, input_files, metadata, output_files):
        """
        The main function to run makeDesignFiles.

        Parameters:
        ----------

        input_files: dict
            designDir : path to the designDir containin .rmap and .baitmap files
        metadata: dict
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

        commands_params = self.get_design_params(self.configuration)

        logger.info("makeDesignFiles command parameters " + " ".join(commands_params))

        results = self.makeDesignFiles(input_files["RMAP"],
                                       input_files["BAITMAP"],
                                       output_files[".nbpb"],
                                       output_files[".npb"],
                                       output_files[".poe"],
                                       commands_params)

        results = compss_wait_on(results)

        out_design_dir = "".join(input_files["RMAP"].split("/")[-1])
        out_design_dir = "".join(out_design_dir.split(".")[0])

        output_metadata = {
            ".nbpb" : Metadata(
                data_type=".nbpb",
                file_type=".nbpb",
                file_path=out_design_dir+".nbpb",
                sources=[
                    metadata["RMAP"].file_path,
                    metadata["BAITMAP"].file_path
                    ],
                taxon_id=[
                    metadata["RMAP"].taxon_id,
                    metadata["BAITMAP"].taxon_id
                    ],
                meta_data={
                    "tool" : "makeDesignFiles.py"
                }
            ),
            ".npb" : Metadata(
                data_type=".npb",
                file_type=".npb",
                file_path=out_design_dir+".npb",
                sources=[
                    metadata["RMAP"].file_path,
                    metadata["BAITMAP"].file_path
                    ],
                taxon_id=[
                    metadata["RMAP"].taxon_id,
                    metadata["BAITMAP"].taxon_id
                    ],
                meta_data={
                    "tool" : "makeDesignFiles.py"
                }
            ),
            ".poe" : Metadata(
                data_type=".poe",
                file_type=".poe",
                file_path=out_design_dir+".poe",
                sources=[
                    metadata["RMAP"].file_path,
                    metadata["BAITMAP"].file_path
                    ],
                taxon_id=[
                    metadata["RMAP"].taxon_id,
                    metadata["BAITMAP"].taxon_id
                    ],
                meta_data={
                    "tool" : "makeDesignFiles.py"
                }
            ),
        }

        return output_files, output_metadata
