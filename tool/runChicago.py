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
import shlex
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
    def chicago(self, input_files, output_files, params):
        
        """
        Run and annotate the Capture-HiC peaks. Chicago will create 4 folders under the outpu_prefix 
        folder:
            data :
                output_index.Rds : chicago dara saved on Rds format
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

        #output should be separated into prefix and path
        output_prefix = output_files.split("/")[-1]
        output_dir = "/".join(output_files.split("/")[:-1])


        #I have runChicago.R added to PATH in bin so no need to call Rscript
        args = ["runChicago.R", input_files, output_prefix, "--output-dir", output_dir]
    
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

    @staticmethod # is there any reason for this to be static?
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

        print(command_params)
        return command_params
       
    def run(self, input_files, output_files):
        """
        The main function to run chicago for peak calling. The input files
        are .chinput and are transformed from BAM files using bam2chicago.sh
        input files could be just one file or a comma separated files from more
        than one biological replicate. Technical replicates should be pooled to one
        .chinput

        Parameters
        ----------
        input_files : dict
            list of 
        metadata : dict
        output_files: dict with the output path

        Returns
        -------
        output_files : Dict
            List of locations for the output files, 
        output_metadata : Dict
            List of matching metadata dict objects
        """
        
        command_params = self.get_chicago_params(self.configuration)
        logger.info("Chicago command parameters "+ " ".join(command_params))

        results = self.chicago(input_files["input_files"], output_files["output_files"], command_params)

        results = compss_wait_on(results)
        
        """
        output_metadata = { 
            "output" : Metadata(
                data_type="data_type",
                file_type=input_files["export_format"],
                file_path=output_files["output"],
                sources=[
                    input_metadata["input_files"].file_path,
                ],
                taxon_id=input_metadata["input_files"].taxon_id,
                meta_data={
                    "tool": "Chicago, capture Capture-HiC algorithm"
                }
            )
                          }
        """
        return(results)


config = { "chicago_setting_file": "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/PCHiCdata/inst/extdata/sGM12878Settings/sGM12878.settingsFile",
               "chicago_desing_dir": "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/PCHiCdata/inst/extdata/hg19TestDesign",
              # "chicago_print_memory": None,
               "chicago_cutoff": "5",
               "chicago_export_format": "washU_text",
               #"chicago_export_order": None,
              # "chicago_rda": None,
               #"chicago_save_df_only": None,
    "chicago_examples_prox_dist": "1e6" ,
  #"chicago_examples_full_range": None,
   #"chicago_en_feat_files": None,
   #"chicago_en_feat_folder": None,
   "chicago_en_min_dist": "0",
   "chicago_en_max_dist": "1e6",
   #"chicago_en_full_cis_range": None,
   "chicago_en_sample_no": "100",}
   #"chicago_en_trans": None,
    #   "chicago_features_only":None}

output_ = {"output_files": "entuano"}

Metadata = {"input_files": "/Users/pacera/MuG/output"}

input_files = {"input_files": "/Users/pacera/MuG/chicagoTeam-chicago-ceffddda8ea3/PCHiCdata/inst/extdata/GMchinputFiles/GM_rep1.chinput"}

test1 = ChicagoTool(config)


print(test1.run(input_files, output_))

        




"""
input_files, output_prefix, setting_file=None, desing_dir="",
                print_memory=None, cutoff="5", export_format="washU_text", export_order=None, 
                rda=None, save_df_only=None, examples_prox_dist="1e6", examples_full_range=None, 
                output_dir=None, en_feat_files=None, en_feat_list=None, 
                en_feat_folder=None, en_min_dist="0", en_max_dist="1e6", en_full_cis_range=None,
                en_sample_no="100", en_trans=None, features_only=None)



input_files : str
        output_prefix : str
        setting_file : str 
        desing_dir : str
        print_memory : flag
        cutoff : int
        export_format : str
        export_order : str 
            "seqMonk", "interBed", "washU_text", "washU_track" (one or more comma separated)
        rda: flag
        save_df_only : flag
        examples_prox_dist : str
        examples_full_range : flag
        output_dir : str
        en_feat_files : str
            comma separated files with genomic features
        en_feat_list : str
        en_feat_folder: str
        en_min_dist : str
        en_max_dist : str
        en_full_cis_range : flag
        en_sample_no : str
        en_trans: flag 
        features_only: flag

"""


















