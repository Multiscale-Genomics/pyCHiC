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
import os
import argparse
import sys
from basic_modules.metadata import Metadata
from process_CHiC import process_CHiC

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
                                     """
                                     This script runs the CHi-C pipeline using pyCHiC
                                     to normalize and call chromatin interactions.
                                     """)

    OPTIONAL = parser._action_groups.pop()
    REQUIRED = parser.add_argument_group('required arguments')

    #Inputs
    REQUIRED.add_argument("--genome_fa",
                          help=" Reference genome in FASTA format",
                          required=True
                         )
    REQUIRED.add_argument("--probes_fa",
                          help=" FASTA format containing the DNA probes"
                               " used to capture the baits",
                          required=True)
    REQUIRED.add_argument("--fastq1", help=" FASTQ forward reads",
                          required=True)
    REQUIRED.add_argument("--fastq2", help=" FASTQ reverse reads",
                          required=True)
    REQUIRED.add_argument("--bowtie_idx",
                          help=" TAR folder containing all files from a"
                               " bowtie2 indexing of the reference genome",
                          required=True)

    #Arguments
    REQUIRED.add_argument("--RE_name",
                          help="Restriction enzyme name used to perform the "
                               "experiment (case sensitive). Example: HindIII",
                          required=True)

    REQUIRED.add_argument("--RE_sequence",
                          help="Restriction enzyme target sequence, "
                               "using a pipe '^' indicating the cutting place."
                               " Example A^AGCTT",
                          required=True)

    OPTIONAL.add_argument("--execution",
                          help=" execution directory",
                          default=".")


    OPTIONAL.add_argument("--minFragLen",
                          help=" The minimum fragment length. If there is any "
                               "fragment shorter than this length, this would be eliminated.",
                          default="150")

    OPTIONAL.add_argument("--maxFragLen",
                          help=" The maximum fragment length. If there is any "
                               "fragment larger than this length, this would be eliminated.",
                          default="40000")

    OPTIONAL.add_argument("--maxLBrownEst",
                          help=" Maximum distance to calculate the Brownian noise",
                          default="1500000")

    OPTIONAL.add_argument("--binsize",
                          help="The size of the bins (in bps) for pooling restriction fragments. ",
                          default="20000")

    OPTIONAL.add_argument("--hicup_longest",
                          help=" Parameter passed to HiCUP. "
                               "Maximum size from the start of "
                               "the read to the restriction size target. ",
                          default="800")

    OPTIONAL.add_argument("--hicup_shortest",
                          help=" Parameter passed to HiCUP. "
                               "Minimum size from the start of the "
                               "read to the restriction size target. ",
                          default="150")

    OPTIONAL.add_argument("--cutoff",
                          help=" significant cutoff for chromatin interactions",
                          default="5")

    OPTIONAL.add_argument("--export_format",
                          help=" pyCHiC export format: seqMonk, interBed, washU_text",
                          default="washU_text")

    OPTIONAL.add_argument("--order",
                          help=" Order of the output interactions: 'position' or 'score'",
                          default="score")

    OPTIONAL.add_argument("--minNPerBait",
                          help=" Minimum number of interactions per bait",
                          default="250")

    OPTIONAL.add_argument("--filterTopPercent",
                          help="Top percentile of trans interactions that are going to be removed"\
                               "durint the estimation of 'other end' bias factor",
                          default="0.01")

    OPTIONAL.add_argument("--minProxOEPerBin",
                          help="Minimum number of trans interactions for every bin during "\
                               "normalisation of other ends",
                          default="150")

    OPTIONAL.add_argument("--minProxB2BPerBin",
                          help="Minimum number of trans interactions for every bin during "\
                               "normalisation of bait to bait in other ends",
                          default="15")

    OPTIONAL.add_argument("--techNoise_minBaitsPerBin",
                          help="Minimum number of trans interactions for every "\
                               "bin during normalisation of baits",
                          default="150")

    OPTIONAL.add_argument("--brownianNoise_samples",
                          help=" Number of times to estimate the brownian noise",
                          default="1")

    OPTIONAL.add_argument("--brownianNoise_subset",
                          help=" Number of samples to take for the estimation "\
                               "of the brownian noise",
                          default="500")

    OPTIONAL.add_argument("--brownianNoise_seed",
                          help=" Seed to control the samples subseting for the "\
                               "brownian estimation",
                          default=None
                         )

    OPTIONAL.add_argument("--weightAlpha",
                          help=" Alpha parameter for the p value weighting",
                          default="34.1157346557331")

    OPTIONAL.add_argument("--weightBeta",
                          help=" Beta parameter for the p value weighting",
                          default="-2.58688050486759")

    OPTIONAL.add_argument("--weightGamma",
                          help=" Gamma parameter for the p value weighting",
                          default="-17.1347845819659")

    OPTIONAL.add_argument("--weightDelta",
                          help=" Delta parameter for the p value weighting",
                          default="-7.07609245521541")

    OPTIONAL.add_argument("--plot_baits",
                          help=" coma separated list of features from baits to plot",
                          default="None")

    parser._action_groups.append(OPTIONAL) # pylint: disable=protected-access


    ARGS = parser.parse_args()

    GENOME_FA = ARGS.genome_fa
    PROBES_FA = ARGS.probes_fa
    FASTQ1 = ARGS.fastq1
    FASTQ2 = ARGS.fastq2
    BOWTIE = ARGS.bowtie_idx
    EXECUTION = ARGS.execution
    RE_NAME = ARGS.RE_name
    RE_SEQUENCE = ARGS.RE_sequence
    MINFRAGLEN = ARGS.minFragLen
    MAXFRAGLEN = ARGS.maxFragLen
    MAXLBROWNEST = ARGS.maxLBrownEst
    BINSIZE = ARGS.binsize
    HICUP_LONGEST = ARGS.hicup_longest
    HICUP_SHORTEST = ARGS.hicup_shortest
    CUTOFF = ARGS.cutoff
    EXPORT_FORMAT = ARGS.export_format
    ORDER = ARGS.order
    MINNPERBAIT = ARGS.minNPerBait
    FILTERTOPPERCENT = ARGS.filterTopPercent
    MINPROXOEPERBIN = ARGS.minProxOEPerBin
    MINPROXB2BPERBIN = ARGS.minProxB2BPerBin
    TECHNOISE_MINBAITSPERBIN = ARGS.techNoise_minBaitsPerBin
    BROWNIANNOISE_SAMPLES = ARGS.brownianNoise_samples
    BROWNIANNOISE_SUBSET = ARGS.brownianNoise_subset
    BROWNIANNOISE_SEED = ARGS.brownianNoise_seed
    WEIGHTALPHA = ARGS.weightAlpha
    WEIGHTBETA = ARGS.weightBeta
    WEIGHTGAMMA = ARGS.weightGamma
    WEIGHTDELTA = ARGS.weightDelta
    PLOTBAITS = ARGS.plot_baits

    RE_SEQUENCE = RE_SEQUENCE.replace("^", "|")

    INPUT_FILES = {
        "genome_fa" : GENOME_FA,
        "probes_fa" : PROBES_FA,
        "fastq1" : FASTQ1,
        "fastq2" : FASTQ2,
        "bowtie_gen_idx" : BOWTIE
    }

    CONFIGURATION = {
        "execution" : EXECUTION,
        "chic_RE_name" : RE_NAME,
        "chic_RE_sequence" : RE_SEQUENCE,
        "makeDesignFiles_minFragLen" : MINFRAGLEN,
        "makeDesignFiles_maxFragLen" : MAXFRAGLEN,
        "makeDesignFiles_maxLBrownEst" : MAXLBROWNEST,
        "makeDesignFiles_binsize" : BINSIZE,
        "hicup_longest" : HICUP_LONGEST,
        "hicup_shortest": HICUP_SHORTEST,
        "pychic_cutoff" : CUTOFF,
        "pychic_export_format": EXPORT_FORMAT,
        "pychic_order" : ORDER,
        "pychic_minNPerBait" : MINNPERBAIT,
        "pychic_tlb_filterTopPercent" : FILTERTOPPERCENT,
        "pychic_tlb_minProxOEPerBin" : MINPROXOEPERBIN,
        "pychic_tlb_minProxB2BPerBin" : MINPROXB2BPERBIN,
        "pychic_techNoise_minBaitsPerBin": TECHNOISE_MINBAITSPERBIN,
        "pychic_brownianNoise_samples" : BROWNIANNOISE_SAMPLES,
        "pychic_brownianNoise_subset" : BROWNIANNOISE_SUBSET,
        "pychic_brownianNoise_seed": BROWNIANNOISE_SEED,
        "pychic_weightAlpha" : WEIGHTALPHA,
        "pychic_weightBeta" : WEIGHTBETA,
        "pychic_weightGamma" : WEIGHTGAMMA,
        "pychic_weightDelta" : WEIGHTDELTA,
        "pychic_features_plot": PLOTBAITS

    }

    OUTPUT_FILES = {
        "hicup_outdir_tar": "tests/data/test_hicup/output.tar",
        "chinput": "tests/data/test_baitmap/output_chinput.chinput",
        "washU_text" :  "tests/data/test_baitmap/washu_test.txt",
        "pdf_examples": "tests/data/test_baitmap/pdf_examples.pdf",
        "params_out": "tests/data/parameters.txt"
    }

    INPUT_METADATA = {
        "fastq1": Metadata(
            data_type="text",
            file_type="fastq",
            file_path=INPUT_FILES["fastq1"],
            sources="",
            taxon_id=9606,
            meta_data=""
            ),
        "fastq2": Metadata(
            data_type="text",
            file_type="fastq",
            file_path=INPUT_FILES["fastq2"],
            sources="",
            taxon_id=9606,
            meta_data=""
            ),
        "genome_fa" : Metadata(
            data_type="text",
            file_type="fasta",
            file_path=INPUT_FILES["genome_fa"],
            sources="",
            taxon_id=9606,
            meta_data={
                "assembly" : "GRCh38"
                },
            ),
        "bowtie_gen_idx" : Metadata(
            "index_bwa", "", INPUT_FILES["genome_fa"],
            {
                "assembly": "test",
                "tool": "bwa_indexer"
            }
            ),
        "probes_fa" : Metadata(
            "C-HiC probes", "fasta", "tests/data/test_baitmap/baits.fa",
            None, None, 9606),

        "hicup_outdir_tar" : Metadata(
            "TAR", "CHiC_data", "tests/data/SRR3535023_1_2.hicup.bam",
            {"fastq1" : "SRR3535023_1.fastq",
             "fastq2" : "SRR3535023_2.fastq", "genome" : "human_hg19"},
            9606),
        "chinput" : Metadata(
            "data_chicago", "chinput", [], None, None, 9606)
        }

    sys._run_from_cmdl = True # pylint: disable=protected-access

    if not os.path.isdir(CONFIGURATION["execution"]):
        os.makedirs(CONFIGURATION["execution"])


    CHIC_HDL = process_CHiC(CONFIGURATION)
    CHIC_HDL.run(INPUT_FILES, INPUT_METADATA, OUTPUT_FILES)
