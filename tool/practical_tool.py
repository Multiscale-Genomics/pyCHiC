from __future__ import print_function

import os
import shlex
import shutil
import subprocess
import sys
import tarfile

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

 from utils.dummy_pycompss import FILE_IN, FILE_OUT # pylint: disable=ungrouped-imports
 from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
 from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool
from basic_modules.metadata import Metadata

# ------------------------------------------------------------------------------

class bwaIndexerTool(Tool):
    """
    Tool for running indexers over a genome FASTA file
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("BWA Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def bwa_index_genome(self, genome_file):
        """
        Create an index of the genome FASTA file with BWA. These are saved
        alongside the assembly file. If the index has already been generated
        then the locations of the files are returned

        Parameters
        ----------
        genome_file : str
            Location of the assembly file in the file system

        Returns
        -------
        amb_file : str
            Location of the amb file
        ann_file : str
            Location of the ann file
        bwt_file : str
            Location of the bwt file
        pac_file : str
            Location of the pac file
        sa_file : str
            Location of the sa file

        """
        command_line = 'bwa index ' + genome_file

        amb_name = genome_file + '.amb'
        ann_name = genome_file + '.ann'
        bwt_name = genome_file + '.bwt'
        pac_name = genome_file + '.pac'
        sa_name = genome_file + '.sa'

        if os.path.isfile(bwt_name) is False:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()

        return (amb_name, ann_name, bwt_name, pac_name, sa_name)

    @task(file_loc=FILE_IN, idx_out=FILE_OUT)
    def bwa_indexer(self, file_loc, idx_out): # pylint: disable=unused-argument
        """
        BWA Indexer

        Parameters
        ----------
        file_loc : str
            Location of the genome assembly FASTA file
        idx_out : str
            Location of the output index file

        Returns
        -------
        bool
        """
        amb_loc, ann_loc, bwt_loc, pac_loc, sa_loc = self.bwa_index_genome(file_loc)

        # tar.gz the index
        print("BS - idx_out", idx_out, idx_out.replace('.tar.gz', ''))
        idx_out_pregz = idx_out.replace('.tar.gz', '.tar')

        index_dir = idx_out.replace('.tar.gz', '')
        os.mkdir(index_dir)

        idx_split = index_dir.split("/")

        shutil.move(amb_loc, index_dir)
        shutil.move(ann_loc, index_dir)
        shutil.move(bwt_loc, index_dir)
        shutil.move(pac_loc, index_dir)
        shutil.move(sa_loc, index_dir)

        index_folder = idx_split[-1]

        tar = tarfile.open(idx_out_pregz, "w")
        tar.add(index_dir, arcname=index_folder)
        tar.close()

        command_line = 'pigz ' + idx_out_pregz
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        return True

    def run(self, input_files, metadata, output_files):
        """
        Function to run the BWA over a genome assembly FASTA file to generate
        the matching index for use with the aligner

        Parameters
        ----------
        input_files : list
            List containing the location of the genome assembly FASTA file
        meta_data : list
        output_files : list
            List of outpout files generated

        Returns
        -------
        output_files : dict
            index : str
                Location of the index file defined in the input parameters
        output_metadata : dict
            index : Metadata
                Metadata relating to the index file
        """
        results = self.bwa_indexer(
            input_files["genome"],
            output_files["index"]
        )
        results = compss_wait_on(results)

        output_metadata = {
            "index": Metadata(
                data_type="sequence_mapping_index_bwa",
                file_type="TAR",
                file_path=output_files["index"],
                sources=[metadata["genome"].file_path],
                taxon_id=metadata["genome"].taxon_id,
                meta_data={
                    "assembly": metadata["genome"].meta_data["assembly"],
                    "tool": "bwa_indexer"
                }
            )
        }

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------

