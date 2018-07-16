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
import tarfile
from shutil import rmtree
from shutil import move


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

        Are you developing 3D/4D genomic analysis tools?
        join us in Septemeber and learn how to integrate it into a Virtual research
        enviroment.
        """

        print("Running Chicago")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def untar_chinput(chinput_tar):
        """
        This function take as input the tar chinput

        Parameters
        ----------
        chinput_tar: str
            path to the tar file, the tar files should have the same prefix name as the
            tar file

        Returns
        -------
        list of untar files
        """
        if chinput_tar.split(".")[-1] == "tar":
            tar = tarfile.open("chinput_tar")
            tar.extractall()
            tar.close()

            directory_path = os.path.split(chinput_tar)

            if directory_path[0] == "":
                files_dir = os.listdir(".")
            else:
                files_dir = os.listdir(directory_path[0])

            files_dir = [fl for fl in files_dir if fl.split(".")[-1] == ".chinput"]

            return files_dir

        else:
            return chinput_tar

    @task(returns=bool, input_files=FILE_IN, output_prefix=IN, output=FILE_OUT,
          params=IN, RMAP=FILE_IN, BAITMAP=FILE_IN, nbpb=FILE_IN, npb=FILE_IN,
          poe=FILE_IN, setting_file=FILE_IN)
    def chicago(self, input_files, output_prefix, output, params, RMAP,
               BAITMAP, nbpb, npb, poe, setting_file):
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

        Parameters
        ----------
        input_files: str ot comma separated list if there is more than one replicate
        output_prefix: str
        output_dir: str (whole path for the output)
        params: dict

        Returns
        -------
        bool
            writes the output files in the defined location
        """
        output_dir = os.path.split(output)[0]

        script = os.path.join(os.path.dirname(__file__), "scripts/runChicago.R")

        input_files = self.untar_chinput(input_files)

        #check if there are more than one .chinput files
        if isinstance(input_files, list):
            args = ["Rscript", script, ", ".join(input_files),
                    output_prefix,
                    "--output-dir", output_dir,
                    "--settings-file", setting_file]

            args += params

        #I have runChicago.R added to PATH in bin so no need to call Rscript
        else:
            args = ["/usr/local/bin/Rscript", script,
                    input_files,
                    output_prefix,
                    "--output-dir", output_dir,
                    "--settings-file", setting_file]

            args += params

        logger.info("chicago CMD: " + " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:
            tar = tarfile.open(os.path.split(output)[1], "w")
            tar.add(output_dir+"/data",
                    arcname="data")

            tar.add(output_dir+"/diag_plots",
                    arcname="diag_plots")

            tar.add(output_dir+"/examples",
                    arcname="examples")

            tar.add(output_dir+"/enrichment_data",
                    arcname="enrichment_data")
            tar.close()

            rmtree(output_dir+"/data")
            rmtree(output_dir+"/diag_plots")
            rmtree(output_dir+"/examples")
            rmtree(output_dir+"/enrichment_data")

            move(os.path.split(output)[1], output)

            logger.info("Tar folder with chinput output file")


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
            "chicago_design_dir": ["--design-dir", True],
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
            if param in command_parameters and params[param] != "None":
                if command_parameters[param][1]:
                    command_params += [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params += [command_parameters[param][0]]

        return command_params

    def run(self, input_files, metadata, output_files):
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
        output_dir = os.path.split(output_files["output"])[0]
        print("output", output_dir)
        if not os.path.exists(output_dir):
            logger.info("creating output directory: "+output_dir)
            os.makedirs(output_dir)


        command_params = self.get_chicago_params(self.configuration)

        logger.info("Chicago command parameters "+ " ".join(command_params))

        design_dir = os.listdir(self.configuration["chicago_design_dir"])
        print(design_dir)

        RMAP = "".join([fl for fl in design_dir if fl.split(".")[-1] == "rmap"])
        BAITMAP = "".join([fl for fl in design_dir if fl.split(".")[-1] == "baitmap"])
        nbpb = "".join([fl for fl in design_dir if fl.split(".")[-1] == "nbpb"])
        npb = "".join([fl for fl in design_dir if fl.split(".")[-1] == "npb"])
        poe = "".join([fl for fl in design_dir if fl.split(".")[-1] == "poe"])

        results = self.chicago(input_files["chinput"],
                               self.configuration["chicago_out_prefix"],
                               output_files["output"],
                               command_params,
                               self.configuration["chicago_design_dir"]+"/"+RMAP,
                               self.configuration["chicago_design_dir"]+"/"+BAITMAP,
                               self.configuration["chicago_design_dir"]+"/"+nbpb,
                               self.configuration["chicago_design_dir"]+"/"+npb,
                               self.configuration["chicago_design_dir"]+"/"+poe,
                               input_files["setting_file"])

        results = compss_wait_on(results)

        output_metadata = {
            "output" : Metadata(
                data_type="tar",
                file_type="tar",
                file_path=output_files["output"],
                sources=[
                    metadata["chinput"].file_path,
                ],
                taxon_id=metadata["chinput"].taxon_id,
                meta_data={
                    "tool": "run_chicago"
                }
            )
        }

        return output_files, output_metadata
