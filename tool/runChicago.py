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
    def chicago(self, input_files, output_prefix, setting_file=None, desing_dir="",
                print_memory=None, cutoff="5", export_format="washU_text", export_order=None, 
                rda=None, save_df_only=None, examples_prox_dist="1e6", examples_full_range=None, 
                output_dir=None, en_feat_files=None, en_feat_list=None, 
                en_feat_folder=None, en_min_dist="0", en_max_dist="1e6", en_full_cis_range=None,
                en_sample_no="100", en_trans=None, features_only=None):
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
        examples_prox_dist : int
        examples_full_range : flag
        output_dir : str
        en_feat_files : str
            comma separated files with genomic features
        en_feat_list : str
        en_feat_folder: str
        en_min_dist : str
        en_max_dist : int
        en_full_cis_range : flag
        en_sample_no : int
        en_trans: flag 
        features_only: flag

        Returns:
        --------
        bool 
            writes the output files in the defined location
        
        """
        #select variables that are not None
        total_args = [arg for arg in locals()]
        print(total_args)
        selected_var = [arg for arg in total_args if locals()[arg] != None]
        
        flags = ["print_memory",
                 "save_df_only", 
                 "examples_full_range",
                 "en_full_cis_range",
                 "en_trans",
                 "features_only"
                 ]

        #dictionary containin the flags from the arguments
        arg_flags = {"setting_file": "--settings-file",
                     "desing_dir": "--design-dir",
                     "print_memory": "--print-memory",
                     "cutoff": "--cutoff",
                     "export_format": "--export-format",
                     "export_order": "--export-order",
                     "rda": "--rda",
                     "save_df_only": "--save-df-only",
                     "examples_prox_dist": "--examples-prox-dist",
                     "examples_full_range": "--examples-full-range",
                     "output_dir": "--output-dir",
                     "en_feat_files": "--en-feat-files",
                     "en_feat_list": "--en-feat-list",
                     "en_feat_folder" :"--en-feat-folder",
                     "en_min_dist": "--en-min-dist",
                     "en_max_dist": "--en-max-dist",
                     "en_full_cis_range": "--en-full-cis-range",
                     "en_sample_no": "--en-sample-no",
                     "en_trans": "--en-trans",
                     "features_only": "--features-only"
                     }
        
        rscript = os.path.join(os.path.dirname(__file__), "../scripts/runChicago.R")

        #include parametres that are not None to call the Rscript
        args = ["Rscript", rscript, input_files, output_prefix]
        for arg in selected_var:
            if arg not in ["self", "input_files", "output_prefix"]:
                if arg in flags:
                    args.append(arg_flags[arg])
                else:
                    args.append(arg_flags[arg])
                    args.append(locals()[arg])
       
        print(args)
        logger.info("chicago CMD: " + " ".join(arg))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            if locals()[output_dir] != None: 
                with open(output_dir + "/data/"+output_prefix+".Rds", "r") as file_in:
                    pass
            else:
                with open(output_prefix + "/data/"+output_prefix+".Rds", "r") as file_in:
                    pass
        except IOError:
            logger.fatal("chicago failed to generate peak file")
            logger.fatal("chicago stdout" + proc_out)
            logger.fatal("chicago stderr" + proc_err)
            return False

        return True


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
            list of 
        metadata : dict

        Returns
        -------
        output_files : Dict
            List of locations for the output files.
        output_metadata : Dict
            List of matching metadata dict objects
        """
        
        results = self.chicago(input_files["input_files"],
                               input_files["output_prefix"],
                               input_files["setting_file"],
                               input_files["desing_dir"],
                               input_files["print_memory"],
                               input_files["cutoff"],
                               input_files["export_format"],
                               input_files["export_order"],
                               input_files["rda"],
                               input_files["save_df_only"],
                               input_files["examples_prox_dist"],
                               input_files["examples_full_range"],
                               output_files["output_dir"],
                               input_files["en_feat_files"],
                               input_files["en_feat_files"],
                               input_files["en_feat_folder"],
                               input_files["en_min_dist"],
                               input_files["en_max_dist"],
                               input_files["en_full_cis_range"],
                               input_files["en_sample_no"],
                               input_files["en_trans"],
                               input_files["features_only"])

        results = compss_wait_on(results)
        
        output_metadata = { 
    
            "output" : Metadata(
                data_type="data_type",
                file_type=input_files["export_format"],
                file_path=output_files["output_dir"],
                sources=[
                    input_metadata["input_files"].file_path,
                ],
                taxon_id=input_metadata["input_files"].taxon_id,
                meta_data={
                    "tool": "Chicago, capture Capture-HiC algorithm"
                }
            )
                          }

        return(results, output_metadata)
























