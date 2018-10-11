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
from basic_modules.metadata import Metadata

from process_CHiC import process_CHiC # pylint: disable=import-error

def test_process_CHiC():
    """
    Test for the process_chicago_CHiC
    """

    path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "fastq1" : path + "test_truncater/SRR3535023_1.fastq",
        "fastq2" : path + "test_truncater/SRR3535023_2.fastq",
        "genome_fa" : path + "test_baitmap/chr21_hg19.fa",
        "genome_idx" : path + "test_baitmap/bwa.tar.gz",
        "probes_fa" : path + "test_baitmap/baits.fa",
        "Rtree_file_dat" : path + "test_rmap/rtree_file.dat",
        "Rtree_file_idx" : path + "test_rmap/rtree_file.idx",
        "chr_handler" : path + "test_baitmap/chr_handler.txt",
        "bowtie_gen_idx" : path + "test_baitmap/chr21_hg19.fa.bt2.tar.gz",
        "RMAP" : path + "test_run_chicago/test.rmap",
        "BAITMAP": path + "test_run_chicago/test.baitmap",
        "hicup_outdir_tar" : path + "test_hicup/output.tar",
        "chinput": [
            path + "test_run_chicago/data_chicago/GM_rep1.chinput",
            path + "test_run_chicago/data_chicago/GM_rep2.chinput"
            ],
        "setting_file" : path + "test_run_chicago/data_chicago/sGM12878.settingsFile",
        "rmap_chicago" : path + "test_run_chicago/data_chicago/h19_chr20and21.rmap",
        "baitmap_chicago" : path + "test_run_chicago/data_chicago/h19_chr20and21.baitmap",
        "nbpb_chicago" : path + "test_run_chicago/data_chicago/h19_chr20and21.nbpb",
        "poe_chicago" : path + "test_run_chicago/data_chicago/h19_chr20and21.poe",
        "npb_chicago" : path + "test_run_chicago/data_chicago/h19_chr20and21.npb",

    }


    output_files = {
        "RMAP" : path + "test_run_chicago/test.rmap",
        "Rtree_file_dat" : path + "test_rmap/rtree_file.dat",
        "Rtree_file_idx" : path + "test_rmap/rtree_file.idx",
        "bait_sam" :  path + "test_baitmap/baits.sam",
        "out_bam" : path +  "test_baitmap/baits.bam",
        "out_baitmap" : path + "test_run_chicago/test.baitmap",
        "hicup_outdir_tar" : path + "test_hicup/output.tar",
        "nbpb" : path + "test_run_chicago/test.nbpb",
        "npb" : path + "test_run_chicago/test.npb",
        "poe" : path + "test_run_chicago/test.poe",
        "chinput" :  path + "test_bam2chicago_tool/output_chinput.chinput",
        "output": path + "test_run_chicago/data_chicago/out_run_chicago.tar",
        "chr_handler" : path + "test_baitmap/chr_handler.txt"
    }

    configuration = {
        "renzime" : {"HindIII" : 'A|AGCTT'},
        "hicup_renzyme" : "A^AGCTT,HindIII",
        "genome_name" : "test_hg19",
        "hicup_bowtie2_loc": "/usr/bin/bowtie2",
        "hicup_longest": "800",
        "hicup_shortest": "150",
        "hicup_outdir": path + "test_hicup/output",
        "hicup_zip": "True",
        "makeDesignFiles_minFragLen" : "150",
        "makeDesignFiles_maxFragLen" : "40000",
        "makeDesignFiles_maxLBrownEst" : "1500000",
        "makeDesignFiles_binsize" : "20000",
        "makeDesignFiles_removeb2b" : True,
        "makeDesignFiles_removeAdjacent" : True,
        "makeDesignFiles_outfilePrefix" : path + "test_run_chicago/test",
        #"makeDesignFiles_designDir" : path + "test_run_chicago",
        "makeDesignFiles_rmap" : path + "test_run_chicago/test.rmap",
        "makeDesignFiles_baitmap" :  path + "test_run_chicago/test.baitmap",
        "chicago_design_dir": path + "/test_run_chicago/data_chicago",
        "chicago_print_memory": "None",
        "chicago_out_prefix" : "output_test",
        "chicago_cutoff": "5",
        "chicago_export_format": "washU_text",
        "chicago_export_order": "None",
        "chicago_rda": "None",
        "chicago_save_df_only": "None",
        "chicago_examples_prox_dist": "1e6",
        "chicago_examples_full_range": "None",
        "chicago_en_feat_files": "None",
        "chicago_en_min_dist": "0",
        "chicago_en_max_dist": "1e6",
        "chicago_en_full_cis_range": "None",
        "chicago_en_sample_no": "100",
        "chicago_en_trans": "None",
        "chicago_features_only": "None"
        }

    input_metadata = {
        "fastq1": Metadata(
            data_type="text",
            file_type="fastq",
            file_path=input_files["fastq1"],
            sources="",
            taxon_id=9606,

            meta_data=""
            ),
        "fastq2": Metadata(
            data_type="text",
            file_type="fastq",
            file_path=input_files["fastq2"],
            sources="",
            taxon_id=9606,
            meta_data=""
            ),
        "genome_fa" : Metadata(
            data_type="text",
            file_type="fasta",
            file_path=input_files["genome_fa"],
            sources="",
            taxon_id=9606,
            meta_data="GRCh38",
            ),
        "genome_idx" : Metadata(
            "index_bwa", "", input_files["genome_fa"],
            {
                "assembly": "test",
                "tool": "bwa_indexer"
            }
            ),
        "probes_fa" : Metadata(
            "C-HiC probes", "fasta", path + "test_baitmap/baits.fa",
            None, None, 9606),

        "Rtree_file_dat" : Metadata(
            "Rtree files", [".dat", ".idx"], path + "test_rmap/rtree_file",
            {"genome" : path + "test_rmap/chr21_hg19.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606
            ),

        "Rtree_file_idx" : Metadata(

            "Rtree files", [".dat", ".idx"], path + "test_rmap/rtree_file",
            {"genome" : path + "test_rmap/chr21_hg19.fa",
             "RE" : {"HindIII" : 'A|AGCTT'}},
            None, 9606
            ),
        "RMAP" : Metadata(
            "data_chicago_input", ".rmap",
            path + "test_run_chicago", None, {}, 9606),
        "BAITMAP" : Metadata(
            "data_chicago_input", ".baitmap",
            path + "test_run_chicago", None, {}, 9606),

        "hicup_outdir_tar" : Metadata(
            "TAR", "CHiC_data", path + "/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606),
        "chinput" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606)
        }

    chic_hdl = process_CHiC(configuration)
    chic_hdl.run(input_files, input_metadata, output_files)

    #assert makeRmap_Tool.py
    assert os.path.getsize(output_files["Rtree_file_dat"]) > 0
    assert os.path.getsize(output_files["Rtree_file_idx"]) > 0

    #assert makeBaitmap.py
    assert os.path.getsize(output_files["out_bam"]) > 0
    assert os.path.getsize(output_files["out_baitmap"]) > 0

    #assert hicup
    assert os.path.isfile(output_files["hicup_outdir_tar"]) is True
    assert os.path.getsize(output_files["hicup_outdir_tar"]) > 0

    assert os.path.isfile(path + "test_run_chicago/test" + ".nbpb") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".nbpb") > 0

    assert os.path.isfile(path + "test_run_chicago/test" + ".npb") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".npb") > 0

    assert os.path.isfile(path + "test_run_chicago/test" + ".poe") is True
    assert os.path.getsize(path + "test_run_chicago/test" + ".poe") > 0

    assert os.path.isfile(output_files["chinput"]) is True
    assert os.path.getsize(output_files["chinput"]) > 0

    assert os.path.isfile(output_files["output"]) is True

    assert os.path.getsize(output_files["output"]) > 0
