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
import tarfile

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

from tool.aligner_utils import alignerUtils
from tool.common import common

##############################################################


class hicup(Tool):
    """
    Tool to run hicup, from fastq to bam files
    """

    def __init__(self, configuration=None):
        """
        Initialising the function

        Parameters
        ----------
        configuration: dict
        """

        print("hicup initialising")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, genome_file_name=IN, genome_idx=FILE_IN,
          bt2_1_file=FILE_OUT, bt2_2_file=FILE_OUT, bt2_3_file=FILE_OUT,
          bt2_4_file=FILE_OUT, bt2_rev1_file=FILE_OUT, bt2_rev2_file=FILE_OUT)
    def untar_index(  # pylint: disable=too-many-locals,too-many-arguments
            self, genome_file_name, genome_idx,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,
            bt2_rev1_file, bt2_rev2_file):
        """
        Extracts the Bowtie2 index files from the genome index tar file.
        Parameters
        ----------
        genome_file_name : str
            Location string of the genome fasta file
        genome_idx : str
            Location of the Bowtie2 index file
        bt2_1_file : str
            Location of the <genome>.1.bt2 index file
        bt2_2_file : str
            Location of the <genome>.2.bt2 index file
        bt2_3_file : str
            Location of the <genome>.3.bt2 index file
        bt2_4_file : str
            Location of the <genome>.4.bt2 index file
        bt2_rev1_file : str
            Location of the <genome>.rev.1.bt2 index file
        bt2_rev2_file : str
            Location of the <genome>.rev.2.bt2 index file
        Returns
        -------
        bool
            Boolean indicating if the task was successful
        """
        if "no-untar" in self.configuration and self.configuration["no-untar"] is True:
            return True

        gfl = genome_file_name.split("/")
        au_handle = alignerUtils()
        au_handle.bowtie2_untar_index(
            gfl[-1], genome_idx,
            bt2_1_file, bt2_2_file, bt2_3_file, bt2_4_file,
            bt2_rev1_file, bt2_rev2_file)

        return True

    @staticmethod
    def get_hicup_params(params):
        """
        Function to handle to extraction of commandline parameters and formatting
        them for use with hicup

        Parameters
        ----------
        params : dict
            --bowtie        Specify the path to Bowtie
            --bowtie2       Specify the path to Bowtie 2
            --config        Specify the configuration file
            --digest        Specify the digest file listing restriction
                            fragment co-ordinates
            --example       Produce an example configuration file
            --format        Specify FASTQ format
                            Options: Sanger, Solexa_Illumina_1.0,
                            Illumina_1.3, Illumina_1.5
            --help          Print help message and exit
            --index         Path to the relevant reference genome
                            Bowtie/Bowtie2 indices
            --keep          Keep intermediate pipeline files
            --longest       Maximum allowable insert size (bps)
            --nofill        Hi-C protocol did NOT include a fill-in of
                            sticky ends prior to ligation step and
                            therefore FASTQ reads shall be truncated
                            at the Hi-C restriction enzyme cut site
                            (if present) sequence is encountered
            --outdir        Directory to write output files
            --quiet         Suppress progress reports (except
                            warnings)
            --shortest      Minimum allowable insert size (bps)
            --temp          Write intermediate files (i.e. all except
                            summaryfiles and files generated by HiCUP
                            Deduplicator) to a specified directory
            --threads       Specify the number of threads, allowing
                            simultaneous
                            processing of multiple files
            --version       Print the program version and exit
            --zip           Compress output

        Returns
        -------
        list
        """
        command_params = []


        command_parameters = {
            "hicup_bowtie2_loc": ["--bowtie2", True],
            "hicup_bowtie_loc": ["--bowtie", True],
            "hicup_genome_digest" : ["--digest", True],
            "hicup_index" : ["--index", True],
            "hicup_longest": ["--longest", True],
            "hicup_shortest": ["--shortest", True],
            "hicup_outdir": ["--outdir", True],
            "hicup_format": ["--format", False],
            "hicup_keep": ["--keep", False],
            "hicup_nofill": ["--nofill", False],
            "hicup_quite": ["--quiet", False],
            "hicup_temp": ["--temp", False],
            "hicup_version": ["--version", False],
            "hicup_zip": ["--zip", False]
            }

        for param in params:
            if param in command_parameters and params[param] != "None":
                if command_parameters[param][1]:
                    if command_parameters[param][0] == "--outdir":
                        continue
                    if command_parameters[param][0] == "--bowtie2":
                        continue

                    command_params += [command_parameters[param][0], params[param]]
                else:
                    if command_parameters[param][0]:
                        command_params += [command_parameters[param][0]]

        return command_params

    def digest_genome(self, genome_name, re_enzyme, genome_loc, re_enzyme2):
        """
        This function takes a genome and digest it using a restriction enzyme
        specified

        Parameters
        ----------
        genome_name: str
            name of the output genome

        re_enzyme: str
            name of the enzyme used to cut the genome
            format example A^GATCT,BglII .
        genome_loc: str
            location of the genome in FASTA format
        re_enzyme2: str
            Restriction site 2 refers to the second,
            optional (other DNA shearing techniques such as
            sonication may be used) enzymatic digestion.
            This restriction site does NOT form a Hi-C ligation junction.
            This is the restriction enzyme that is used when the Hi-C
            sonication protocol is not followed. Typically the sonication
            protocol is followed.
        """
        if re_enzyme2 == "enzyme2":
            args = ["hicup_digester",
                    "--genome", genome_name,
                    "--re1", re_enzyme,
                    genome_loc]
        else:
            args = ["hicup_digester",
                    "--genome", genome_name,
                    "--re1", re_enzyme,
                    "--re2", re_enzyme2,
                    genome_loc]
        try:
            logger.info("hicup_digester CMD: " + " ".join(args))

            process = subprocess.Popen(
                ' '.join(args),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
                )

            process.wait()

        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, args))

        files_dir = os.listdir(".")

        digest_genome = [file_ for file_ in files_dir if \
            file_.startswith("Digest_"+genome_name)]

        return "".join(digest_genome)


    def hicup_alig_filt(self, params, genome_digest, genome_index,
                        genome_loc, fastq1, fastq2, outdir_tar, input_files):
        """
        This function aling the HiC read into a reference
        genome and filter them

        Parameters
        ----------
        bowtie2_loc:
        genome_index: str
            location of genome indexed with bowtie2
        digest_genome: str
            location of genome digested
        fastq1: str
            location of fastq2 file
        fastq2: str
            location of fastq2

        Returns
        -------
        Bool
        """
        folder = os.path.split(outdir_tar)[0]+"/"+ \
                os.path.split(outdir_tar)[1].split(".")[0]

        if os.path.isdir(folder) is False:
            os.mkdir(folder)

        if "genome_fa_public" in input_files:
            genome_loc = "/orozco/services/MuGdev/MuG_public/refGenomes/hg19/BOWTIE2/bt2/hg19.fa"

        else:
            index_files = {
                "1.bt2": genome_loc + ".1.bt2",
                "2.bt2": genome_loc + ".2.bt2",
                "3.bt2": genome_loc + ".3.bt2",
                "4.bt2": genome_loc + ".4.bt2",
                "rev.1.bt2": genome_loc + ".rev.1.bt2",
                "rev.2.bt2": genome_loc + ".rev.2.bt2"
            }

            logger.progress("Untar Index: "+genome_loc+", "+genome_index)
            self.untar_index(
                genome_loc,
                genome_index,
                index_files["1.bt2"],
                index_files["2.bt2"],
                index_files["3.bt2"],
                index_files["4.bt2"],
                index_files["rev.1.bt2"],
                index_files["rev.2.bt2"]
                )

        hicup_args = [
            "hicup",
            "--index", genome_loc,
            "--digest", genome_digest,
            fastq1,
            fastq2
            ]

        hicup_args = hicup_args + params + ["--bowtie2", "/home/compss/bin/bowtie2",
                                            "--threads", "2",
                                            "--outdir", folder]

        logger.info("arguments for hicup:" + " ".join(hicup_args))


        try:
            process = subprocess.Popen(" ".join(hicup_args), shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()

            logger.info("TARING output folder")

            tar_file = outdir_tar
            archive_name = os.path.split(outdir_tar)[1].split(".")[0]

            onlyfiles = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

            tar = tarfile.open(tar_file, "w")

            for tmp_file in onlyfiles:
                tar.add(
                    os.path.join(folder, tmp_file),
                    arcname=os.path.join(archive_name, tmp_file)
                )

                os.remove(os.path.join(folder, tmp_file))
            tar.close()

            shutil.rmtree(folder)

            shutil.move(tar_file, outdir_tar)

            for indexed_file in index_files:
                os.remove(index_files[indexed_file])

            return True

        except IOError:
            return False

    @task(returns=bool,
          params=IN,
          genome_digest=FILE_IN,
          genome_index=FILE_IN,
          genome_loc=FILE_IN,
          fastq1=FILE_IN,
          fastq2=FILE_IN,
          outdir_tar=FILE_OUT,
          input_files=IN)
    def hicup_alig_filt_runner(self, params, genome_digest, genome_index,
                               genome_loc, fastq1, fastq2, outdir_tar, input_files):
        """
        This function runs the hicup_alig_filt

        Parameters
        ----------
        bowtie2_loc:
        genome_index: str
            location of genome indexed with bowtie2
        digest_genome: str
            location of genome digested
        fastq1: str
            location of fastq2 file
        fastq2: str
            location of fastq2

        Returns
        -------
        Bool
        """
        self.hicup_alig_filt(params, genome_digest, genome_index,
                             genome_loc, fastq1, fastq2, outdir_tar, input_files)

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Function that runs and pass the parameters for
        all the functions

        Parameters
        ----------
        input_files: dict
        metadata: dict
        output_files: dict
        """
        if "genome_fa_public" in input_files:
            input_files["genome_fa"] = input_files.pop("genome_fa_public")
            input_metadata["genome_fa"] = input_metadata.pop("genome_fa_public")

            input_files["genome_idx"] = input_files.pop("genome_idx_public")
            input_metadata["genome_idx"] = input_metadata.pop("genome_idx_public")

            input_files["bowtie_gen_idx"] = input_files.pop("bowtie_gen_idx_public")
            input_metadata["bowtie_gen_idx"] = input_metadata.pop("bowtie_gen_idx_public")

        self.configuration["genome_name"] = "digested_genome"
        self.configuration["hicup_zip"] = "True"

        #check if there is any digest files and delete it
        files_dir = os.listdir(".")
        for file_ in files_dir:
            if file_.startswith("Digest_"+self.configuration["genome_name"]):
                os.remove(file_)


        output_files["hicup_outdir_tar"] = self.configuration["execution"]+"/"+\
                                           os.path.split(output_files["hicup_outdir_tar"])[1]


        self.configuration["hicup_outdir"] = output_files["hicup_outdir_tar"].split(".")[0]


        RE = str(self.configuration["chic_RE_sequence"].replace("|", "^"))
        enzyme = str(self.configuration["chic_RE_name"])

        re_enzyme = RE+","+enzyme

        if "renzyme_name2" in self.configuration:

            genome_d = self.digest_genome(
                self.configuration["genome_name"],
                re_enzyme,
                input_files["genome_fa"],
                self.configuration["renzyme_name2"]
                )
        else:

            genome_d = self.digest_genome(
                self.configuration["genome_name"],
                re_enzyme,
                input_files["genome_fa"],
                "enzyme2"
            )

        parameters_hicup = self.get_hicup_params(self.configuration)

        #if os.path.isdir(self.configuration["hicup_outdir"]) is False:
        #    os.mkdir(self.configuration["hicup_outdir"])

        if "genome_fa_public" in input_files:


        variable = self.hicup_alig_filt(# pylint: disable=too-many-locals,too-many-arguments
            parameters_hicup,
            genome_d,
            input_files["bowtie_gen_idx"],
            input_files["genome_fa"],
            input_files["fastq1"],
            input_files["fastq2"],
            output_files["hicup_outdir_tar"],
            input_files)

        os.remove(genome_d)
        #variable = compss_wait_on(variable)

        output_metadata = {
            "hicup_outdir_tar" : Metadata(
                data_type="data_CHiC",
                file_type="TAR",
                file_path=output_files["hicup_outdir_tar"],
                sources=[
                    input_metadata["genome_fa"].file_path,
                    input_metadata["fastq1"].file_path,
                    input_metadata["fastq1"].file_path
                    ],
                taxon_id=input_metadata["genome_fa"].taxon_id,
                meta_data={"tool": "hicup_tool"}
                )
        }

        return output_files, output_metadata
