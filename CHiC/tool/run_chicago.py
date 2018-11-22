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
import tarfile
from shutil import rmtree
from shutil import move
from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on, compss_delete_file, compss_open
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on, compss_delete_file, compss_open # pylint: disable=ungrouped-imports

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
        ----------
        configuration : dict
            dictionary containing parameters that define how the operation
            should be carried out, which are specific to the tool.

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
            tar = tarfile.open(chinput_tar)
            tar.extractall(path=os.path.split(chinput_tar)[0])
            tar.close()

            directory_path = os.path.split(chinput_tar)[0]+"/chinput"

            files_dir = os.listdir(directory_path)

            chinput_paths = []
            for chinput_name in files_dir:
                chinput_paths.append(directory_path+"/"+chinput_name)

            files_comman_sp = str(",".join(chinput_paths))

            return files_comman_sp

        else:
            return chinput_tar

    @task(returns=bool, input_files=FILE_IN, output=FILE_OUT, params=IN,
          setting_file=FILE_IN, rmap=FILE_IN, baitmap=FILE_IN, nbpb=FILE_IN,
          npb=FILE_IN, poe=FILE_IN, washu=FILE_OUT, pdf=FILE_OUT)
    def chicago(self, input_files, output_prefix, output, params, setting_file,
                rmap, baitmap, nbpb, npb, poe, washu, pdf):
        """
        Run and annotate the Capture-HiC peaks. Chicago will create 4 folders under the outpu_prefix
        data :
        output_index.Rds --> chicago data saved on Rds format
        output_index_params.txt --> parameters used to run Chicago
        output_index.export_format --> chicago output in the chosen format
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

        rlib = os.path.join(os.getcwd(), "tmp_R_lib")
        if not os.path.exists(rlib):
            os.makedirs(rlib)

        input_untared = self.untar_chinput(input_files)

        args = ["Rscript", script,
                input_untared,
                output_prefix,
                "--output-dir", output_dir,
                "--settings-file", setting_file,
                "--design-dir", os.path.split(rmap)[0]]

        args += params
        print(" ".join(args))
        logger.info("chicago CMD: " + " ".join(args))

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        proc_out, proc_err = process.communicate()

        try:

            move(output_dir+"/data/"+output_prefix+"_washU_text.txt",
                 washu)

            move(output_dir+"/examples/"+output_prefix+"_proxExamples.pdf",
                 pdf)

            tar = tarfile.open(output, "w")
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

            logger.info("Tar folder with chinput output file")

            return True

        except IOError:
            logger.fatal("chicago failed to generate peak file")
            logger.fatal("chicago stdout" + proc_out)
            logger.fatal("chicago stderr" + proc_err)

            return False


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
            #"chicago_design_dir": ["--design-dir", True],
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

        print(command_params)
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
        input_metadata : dict
        output_files: dict with the output path

        Returns
        -------
        output_files : Dict
            List of locations for the output files,
        output_metadata : Dict
            List of matching metadata dict objects
        """
        #delete files that are not returned to the user
        rtree_file_dat = "tests/data/test_rmap/rtree_file.dat"
        rtree_file_idx = "tests/data/test_rmap/rtree_file.idx"
        chr_handler = "tests/data/test_baitmap/chr_handler.txt"
        RMAP = "tests/data/test_run_chicago/test.rmap"
        out_baitmap = "tests/data/test_run_chicago/test.baitmap"
        bait_sam = "tests/data/test_baitmap/baits.sam"
        nbpb = "tests/data/test_run_chicago/test.nbpb"
        npb = "tests/data/test_run_chicago/test.npb"
        poe = "tests/data/test_run_chicago/test.poe"
        out_bam = "tests/data/test_baitmap/baits.bam"
        sorted_bam = self.configuration["execution"] + "/" + "sorted_bam"


        self.configuration["chicago_design_dir"] = "tests/data/test_run_chicago/data_chicago"

        #chinput = "tests/data/test_run_chicago/data_chicago/GM_rep1.chinput"
        output_files["chinput"] = self.configuration["execution"]+"/"+\
                                    os.path.split(output_files["chinput"])[1]


        output_files["output"] = self.configuration["execution"]+"/"+\
                                    os.path.split(output_files["output"])[1]

        hicup_folder = self.configuration["execution"]+"/"+\
                                           os.path.split(output_files["hicup_outdir_tar"])[1]


        washu = self.configuration["execution"]+\
            "/"+output_files["washU_text"]

        pdf = self.configuration["execution"]+\
            "/"+output_files["pdf_examples"]



        command_params = self.get_chicago_params(self.configuration)

        logger.info("Chicago command parameters "+ " ".join(command_params))

        """
        washu = self.configuration["execution"]+\
            "/"+self.configuration["chicago_out_prefix"]+"_washU_text.txt"

        pdf = self.configuration["execution"]+\
            "/"+self.configuration["chicago_out_prefix"]+"_proxExamples.pdf"

        """

        results = self.chicago(output_files["chinput"],
                               self.configuration["chicago_out_prefix"],
                               output_files["output"],
                               command_params,
                               input_files["setting_file"],
                               RMAP,
                               out_baitmap,
                               nbpb,
                               npb,
                               poe,
                               washu,
                               pdf
                              )

        compss_delete_file(rtree_file_idx)
        compss_delete_file(rtree_file_dat)
        compss_delete_file(chr_handler)
        compss_delete_file(RMAP)
        compss_delete_file(out_baitmap)
        compss_delete_file(bait_sam)
        compss_delete_file(npb)
        compss_delete_file(nbpb)
        compss_delete_file(poe)
        compss_delete_file(out_bam)
        compss_delete_file(sorted_bam)

        files_dir = os.listdir(self.configuration["execution"])
        for file_ in files_dir:
            if file_.startswith("Digest_"+self.configuration["genome_name"]):
                os.remove(file_)

        output_metadata = {
            "output" : Metadata(
                data_type="chicago_CHIC",
                file_type="TAR",
                file_path=output_files["output"],
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["fastq1"].file_path,
                    input_metadata["fastq2"].file_path
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago",

                }
            ),

            "washU_text" : Metadata(
                data_type="chicago_CHIC",
                file_type="TXT",
                file_path=washu,
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["fastq1"].file_path,
                    input_metadata["fastq2"].file_path
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago",
                }
            ),

            "pdf_examples" : Metadata(
                data_type="chicago_CHIC",
                file_type="PDF",
                file_path=pdf,
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["fastq1"].file_path,
                    input_metadata["fastq2"].file_path
                ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={
                    "tool": "process_CHiC",
                    "tool_description" : "run_chicago",
                }
            )
        }

        return output_files, output_metadata
