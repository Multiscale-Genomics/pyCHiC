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
import pandas as pd
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

###################################################

class ChicagoTool(Tool):
    """
    tool for running the CHiCAGO algorithm
    """

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.

        Parameters
        -----------
        configuration : dict
            dictionary containing parameters that define how the operation
            should be carried out, which are specific to the tool.

        """

        print("Running Chicago")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    #@task(some decorators)
    def chicago(self, input_files, output_prefix, output_dir, params):
        """
        Run and annotate the Capture-HiC peaks. Chicago will create 4 folders under the outpu_prefix
        folder:
            data :
                output_index.Rds : chicago data saved on Rds format
                output_index_params.txt : parameters used to run Chicago
                output_index.export_format : chicago output in the chosen format
            diag_plots :
                3 plots to assest the quality of the output
                (see CHicago Capture-HiC documentation for details)
            enrichment_data:
                files for the feature enrichment output (in case is used)
            examples:
                output_index_proxExamples.pdf: random chosen peaks showing interactions regions
        see http://regulatorygenomicsgroup.org/chicago for more information

        Parameters:
        -----------
        input_files: str ot comma separated list if there is more than one replicate
        output_prefix: str
        output_dir: str (whole path for the output)
        params: dict

        Returns:
        --------
        bool
            writes the output files in the defined location
        """

        #check if there are more than one .chinput files
        if isinstance(input_files, list):
            args = ["../scripts/runChicago.R", ", ".join(input_files),
                    output_prefix, "--output-dir", output_dir]
            args += params

        #I have runChicago.R added to PATH in bin so no need to call Rscript
        else:
            args = ["../scripts/runChicago.R",
                    input_files, output_prefix,
                    "--output-dir", output_dir]

            args += params

        logger.info("chicago CMD: " + " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            os.path.isfile(output_dir+"/data"+output_prefix+".Rds")
        except IOError:
            logger.fatal("chicago failed to generate peak file")
            logger.fatal("chicago stdout" + proc_out)
            logger.fatal("chicago stderr" + proc_err)
            return False

        return True

    @staticmethod
    def get_chicago_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use in the aligner for BWA ALN

        Parameters
        ----------
        params : dict

        Returns
        -------
        list
        """
        command_params = []

        command_parameters = {
            "chicago_setting_file": ["--settings-file", True],
            "chicago_desing_dir": ["--design-dir", True],
            "chicago_print_memory": ["--print-memory", False],
            "chicago_cutoff": ["--cutoff", True],
            "chicago_export_format":["--export-format", True],
            "chicago_export_order": ["--export-order", True],
            "chicago_rda": ["--rda", True],
            "chicago_save_df_only": ["--save-df-only", False],
            "chicago_examples_prox_dist": ["--examples-prox-dist", True],
            "chicago_examples_full_range": ["--examples-full-range", False],
            "chicago_en_feat_files": ["--en-feat-files", True],
            "chicago_en_feat_list": ["--en-feat-list", True],
            "chicago_en_feat_folder" :["--en-feat-folder", True],
            "chicago_en_min_dist": ["--en-min-dist", True],
            "chicago_en_max_dist": ["--en-max-dist", True],
            "chicago_en_full_cis_range": ["--en-full-cis-range", False],
            "chicago_en_sample_no": ["--en-sample-no", True],
            "chicago_en_trans": ["--en-trans", False],
            "chicago_features_only": ["--features-only", False]
            }


        for param in params:
            if param in command_parameters:
                if command_parameters[param][1]:
                    command_params += [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params += [command_parameters[param][0]]

        return command_params

    def run(self, input_files, input_metadata, output_files):
        """
        The main function to run chicago for peak calling. The input files
        are .chinput and are transformed from BAM files using bam2chicago.sh
        input files could be just one file or a comma separated files from more
        than one biological replicate. Technical replicates should be pooled to one
        .chinput

        Parameters
        ----------
        input_files : dict
            list of .chinput files, or str with a single .chinput file
        metadata : dict
        output_files: dict with the output path

        Returns
        -------
        output_files : Dict
            List of locations for the output files,
        output_metadata : Dict
            List of matching metadata dict objects
        """
        #check if the output directory exists, otherwise create it
        if not os.path.exists(output_files["output_dir"]):
            logger.info("creating output directory: "+
                        output_files["output_dir"])
            os.makedirs(output_files["output_dir"])


        command_params = self.get_chicago_params(self.configuration)

        logger.info("Chicago command parameters "+ " ".join(command_params))

        results = self.chicago(input_files["chinput_file"],
                               output_files["output_prefix"],
                               output_files["output_dir"],
                               command_params)

        results = compss_wait_on(results)


        output_metadata = {
            "output" : Metadata(
                data_type="peaks",
                file_type=self.configuration["chicago_export_format"],
                file_path=output_files["output_dir"],
                sources=[
                    input_metadata["chinput_1"].file_path,
                ],
                taxon_id=input_metadata["chinput_1"].taxon_id,
                meta_data={
                    "tool": "Chicago, capture Capture-HiC algorithm"
                }
            )
        }

        return(results, output_metadata)
