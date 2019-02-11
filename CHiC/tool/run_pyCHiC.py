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
from CHiC.tool.pyCHiC import pyCHiC

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
                                     """
                                     This script runs the pyCHiC tool
                                     to normalize and call chromatine interactions.
                                     """)

    OPTIONAL = parser._action_groups.pop()
    REQUIRED = parser.add_argument_group('required arguments')

    #Inputs
    REQUIRED.add_argument("--rmap",
                          help=" rmap file",
                          required=True
                         )
    REQUIRED.add_argument("--baitmap",
                          help=" baitmap file",
                          required=True)
    REQUIRED.add_argument("--nbpb", help=" nbpb file",
                          required=True)
    REQUIRED.add_argument("--npb", help=" npb file",
                          required=True)
    REQUIRED.add_argument("--poe",
                          help=" poe file",
                          required=True)

    REQUIRED.add_argument("--chinput",
                          help="chinput file",
                          required=True)

    OPTIONAL.add_argument("--execution",
                          help="execution directory",
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

    OPTIONAL.add_argument("--cutoff",
                          help=" significant cutoff for chromatin interactions",
                          default="1")

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

    OPTIONAL.add_argument("--output_file",
                          help="output name of the interactions file",
                          default="washU_text.txt")

    OPTIONAL.add_argument("--output_pdf",
                          help="output name for the pdf containing the bait plots",
                          default="pdf_examples.pdf")

    OPTIONAL.add_argument("--output_parameters",
                          help="output name for parameter file",
                          default="parameters.txt")
    OPTIONAL.add_argument("--cpu",
                          default="1")

    parser._action_groups.append(OPTIONAL) # pylint: disable=protected-access


    ARGS = parser.parse_args()

    #inputs
    RMAP = ARGS.rmap
    BAITMAP = ARGS.baitmap
    NBPB = ARGS.nbpb
    NPB = ARGS.npb
    POE = ARGS.poe
    CHINPUT = ARGS.chinput

    #arguments
    EXECUTION = ARGS.execution
    MINFRAGLEN = ARGS.minFragLen
    MAXFRAGLEN = ARGS.maxFragLen
    MAXLBROWNEST = ARGS.maxLBrownEst
    BINSIZE = ARGS.binsize
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
    WASHU = ARGS.output_file
    PDF = ARGS.output_pdf
    PARAMETERS = ARGS.output_parameters
    CPU = ARGS.cpu

    INPUT_FILES = {
        "RMAP" : RMAP,
        "BAITMAP" : BAITMAP,
        "nbpb" : NBPB,
        "npb" : NPB,
        "poe" : POE,
        "chinput" : CHINPUT
    }

    CONFIGURATION = {
        "execution" : EXECUTION,
        "makeDesignFiles_minFragLen" : MINFRAGLEN,
        "makeDesignFiles_maxFragLen" : MAXFRAGLEN,
        "makeDesignFiles_maxLBrownEst" : MAXLBROWNEST,
        "makeDesignFiles_binsize" : BINSIZE,
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
        "pychic_features_plot": PLOTBAITS,
        "pychic_cpu": CPU
    }

    OUTPUT_FILES = {
        "washU_text" :  EXECUTION + "/" + WASHU,
        "pdf_examples": EXECUTION + "/" + PDF,
        "params_out": EXECUTION + "/" + PARAMETERS
    }

    INPUT_METADATA = {
        "chinput": Metadata(
            "data_chicago", "chinput", [], None, None, 9606
            )
    }

    sys._run_from_cmdl = True # pylint: disable=protected-access

    if not os.path.isdir(CONFIGURATION["execution"]):
        os.makedirs(CONFIGURATION["execution"])


    CHIC_HDL = pyCHiC(CONFIGURATION)
    CHIC_HDL.run(INPUT_FILES, INPUT_METADATA, OUTPUT_FILES)
