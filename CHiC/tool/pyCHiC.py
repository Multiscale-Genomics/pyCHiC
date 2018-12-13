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
import sys
import os
import re
import time
import math
from scipy.stats.mstats import gmean
from collections import Counter
from random import sample

import pandas as pd
import numpy as np
from utils import logger
from scipy import stats
from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns

from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2 import robjects

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
from basic_modules.metadata import Metadata # pylint: disable=unused-import

#######################################################

class pyCHiC(Tool): # pylint: disable=invalid-name
    """
    Tool for preprocess the input files
    """

    def __init__(self, configuration=None):
        """
        initialising the function

        Parameters
        ----------
        configuration: dict
         dictionary containing all the arguments and parameters
         to run the tool
        """

        print("CHiCAGO initialising")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @staticmethod
    def exist_file(test_file, file_name):
        """
        Check if the file exists, otherwise send an error message

        Parameter
        ---------
        file: str
            path to the file

        Returns
        -------
        bool
        """
        if isinstance(test_file, list) is True:
            for file in test_file:
                if os.path.isfile(file) is False:
                    logger.fatal(test_file + " does not exists, "
                                 "introduce a valid "+
                                 file_name + " file")
                    return False
        else:
            if os.path.isfile(test_file) is False:
                logger.fatal(test_file + " does not exists, introduce a valid "+
                             file_name + " file")
                return False
        return True

    @staticmethod
    def check_design(design_fl, rmap, baitmap, column_rmap, column_baitmap):
        """
        This function is going to check that the design files
        nbpb, npb and poe have been generated using the given
        rmap and baitmap files

        Parameters
        ----------
        design_fl: str
            path to the design file

        column_rmap: column number of the first line in the
            design file to check for the name of he rmap file
        column_baitmap: column number of the first line in the
            design file to check for the name of he baitmap

        Return
        ------
        Bool
        """
        rmap_name = os.path.split(rmap)[1]
        baitmap_name = os.path.split(baitmap)[1]



        with open(design_fl, "r") as file_in:
            for line in file_in:
                line_hdl = line.rstrip().split("\t")
                if os.path.split(line_hdl[column_rmap].split("=")[1])[1] != rmap_name:
                    logger.fatal(design_fl+" file has not been generated using "+
                                 rmap)
                    logger.fatal("rmap used to generate "+design_fl+" "+
                                 line_hdl[column_rmap].split("=")[1])
                    return False
                if os.path.split(line_hdl[column_baitmap].split("=")[1])[1] != baitmap_name:
                    logger.fatal(design_fl+" file has not been generated using "+
                                 baitmap)
                    logger.fatal("baitmap used to generate "+design_fl+" "+
                                 line_hdl[column_baitmap].split("=")[1])
                    return False
                break
            return True

    def checks(self, pychic_export_format, pychic_export_order, rmap, baitmap,
               nbpb, npb, poe, settingsFile): # pylint: disable=invalid-name
        """
        This fucntion checks that all the input files exists, has the right
        format and have been properly generated.


        Parameters
        ----------
        pychic_export_format: str
            seqMonk", "interBed", "washU_text", "washU_track"
        pychic_export_order: str
            "score" or genomic "position"
        rmap : str
            path to rmap
        baitmap: str
            path to baitmap
        nbpb: str
            path to nbpb
        poe: str
            path to poe
        settingsFile: str
            path to settingsFile

        Returns
        -------
        Bool
        """
        export_formats = ["seqMonk", "interBed", "washU_text", "washU_track"]

        if set(export_formats).intersection(set(pychic_export_format)):
            pass
        else:
            logger.fatal("pychic_export_format must be either seqMonk,"
                         " interBed, washU_text or washU_track"
                         " (or a comma-separated combination of these)")
            return False


        if pychic_export_order not in ["position", "score"]:
            logger.fatal("pychic_export_order must be either"
                         " position (default) or score")
            return False


        logger.info("checking rmap and baitmap files..")
        baitmap_df = pd.read_csv(baitmap, sep="\t", header=None)
        rmap_df = pd.read_csv(rmap, sep="\t", header=None)

        if len(baitmap_df.columns) < 4:
            logger.fatal("There are fewer columns in the baitmapfile than expected."
                         " Check that this file lists the genomic coordinates as well as both "
                         "the IDs and names for each baited fragment")
            logger.fatal("Number of columns found "+str(len(rmap_df.columns)))
            return False

        if len(rmap_df.columns) < 4:
            logger.fatal("There are fewer columns in the rmap file than expected."
                         "This file should have 4 columns, listing the genomic coordinates and "
                         "IDs for each restriction fragment.")
            logger.fatal("Number of columns found "+str(len(rmap_df.columns)))
            return False

        if len(rmap_df.columns) > 4:
            logger.fatal("There are more columns in the rmap file than expected. "
                         "This file should have 4 columns, listing the genomic coordinates "
                         "and IDs for each restriction fragment. Check that rmapfile and baitmap "
                         "files aren't swapped round")
            logger.fatal("Number of columns found "+str(len(rmap_df.columns)))
            return False

        IDbaits = baitmap_df.iloc[:, 3] # pylint: disable=invalid-name
        IDbaitSet = set(IDbaits) # pylint: disable=invalid-name
        if len(IDbaits) != len(IDbaitSet):
            logger.fatal("Duplicated fragment IDs found in baitmapfile. "
                         "Check that the baitmap columns are not swapped "
                         "round (4rd column ID)")
            logger.fatal(str(len(IDbaits)), str(len(IDbaitSet)))

        IDrmap = rmap_df.iloc[:, 3] # pylint: disable=invalid-name
        IDrmapSet = set(IDrmap) # pylint: disable=invalid-name
        if len(IDrmap) != len(IDrmapSet):
            logger.fatal("Duplicated fragment IDs found in baitmapfile. "
                         "Check that the baitmap columns are not swapped "
                         "round (4rd column ID)")
            return False


        if IDbaitSet - IDrmapSet:
            logger.fatal("Error: Some entries of baitmap are not present in "
                         "rmap (listed below).\n")
            for ID in IDbaitSet - IDrmapSet: # pylint: disable=invalid-name
                logger.fatal(ID)
            return False

        logger.info("checking nbpb, npb and poe files...")

        if self.check_design(nbpb, rmap, baitmap, 3, 4) is False:
            return False

        if self.check_design(npb, rmap, baitmap, 7, 8) is False:
            return False

        if self.check_design(poe, rmap, baitmap, 7, 8) is False:
            return False

        return True

    def readSample(self, chinput, bamFile, rmap, baitmap): # pylint: disable=invalid-name
        """
        This function takes the chinput file and filter it according
        with the settings of the experiment

        Parameters
        ----------
        chipnut: str
            path to the chinput file
        bamfile: str
            name of the bamFile used to generate
        rmap: str
            path to the rmap file
        baitmap: str
            path to the baitmap file

        Return
        ------
        chinput_filtered: data frame
        """

        chinput = "".join(chinput)
        logger.info("reading and checking "+chinput)

        rmap_name = os.path.split(rmap)[1]
        baitmap_name = os.path.split(baitmap)[1]

        with open(chinput, "r") as file_in:
            for line in file_in:
                line_hdl = line.rstrip()
                if line_hdl[0] == "#":

                    bam_ = re.compile(r"bamname\S+")
                    rmap_ = re.compile(r"digestfile\S+")
                    baitmap_ = re.compile(r"baitmapfile\S+")

                    bam_chinput = bam_.findall(line_hdl)
                    rmap_chinput = rmap_.findall(line_hdl)
                    baitmap_chinput = baitmap_.findall(line_hdl)

                    if "".join(bam_chinput).split("=")[1] != bamFile:
                        logger.fatal("The chinput file was not generated using "
                                     "the BAM file specified, BAM file provided: "
                                     +bamFile+" BAM file found: "+"".join(bam_chinput).split("=")[1])

                    if "".join(baitmap_chinput).split("=")[1] != baitmap_name:
                        logger.fatal("The chinput file was not generated using "
                                     "the baitmap file specified, baitmap provided: "
                                     +baitmap_name+" baimap found: "+"".join(baitmap_chinput).split("=")[1])

                    if "".join(rmap_chinput).split("=")[1] != rmap_name:
                        logger.fatal("The chinput file was not generated using "
                                     "the rmap file specified, rmap provided"
                                     +rmap_name+" rmap found x"+"".join(rmap_chinput).split("=")[1])
                else:
                    break

        logger.info(chinput+"file checked succesfully")

        x = pd.read_csv(chinput, sep="\t", skiprows=1, dtype= # pylint: disable=invalid-name
                        {"baitID":int,
                         "otherEndID":int,
                         "N":int,
                         "otherEndLen":int,
                         "distSign":float
                        }
                       )

        if len(x.columns) != 5:
            logger.fatal("Wrong number of columns found in chipnut file. "
                         "Chinput file should have 5 columns(baitID, otherEndID,"
                         "N, otherEndLen, distSign)."+ "Number of columns found "+
                         str(x.columns))

        rmap_df = pd.read_csv(rmap, sep="\t", header=None)
        rmap_id = set(rmap_df.iloc[:, 3])

        baitmap_df = pd.read_csv(baitmap, sep="\t", header=None)
        baitmap_id = set(baitmap_df.iloc[:, 3])

        x_rmap = set(x.iloc[:, 1])
        x_baitmap = set(x.iloc[:, 0])

        if x_rmap - rmap_id:
            logger.fatal("Some entries of the chinput file have otherEndID "
                         "(rmap ids) not"
                         "in the rmapfile (listed below). Check that the specified "
                         "design files are correct")
            for ID in x_rmap - rmap_id: # pylint: disable=invalid-name
                logger.fatal(str(ID))

        if x_baitmap - baitmap_id:
            logger.fatal("Some entries of the input file have baitIDs not "
                         "in the baitmap. Check that the specified design files "
                         "are correct")
            for ID in x_baitmap - baitmap_id: # pylint: disable=invalid-name
                logger.info(str(ID))

        logger.info("Filtering interaction with OEs smaller than "+
                    str(self.configuration["pychic_minFragLen"])+" and bigger"
                    "than "+str(self.configuration["pychic_maxFragLen"]))

        x = x.loc[ # pylint: disable=invalid-name
            (x["otherEndLen"] >= self.configuration["pychic_minFragLen"]) &
            (x["otherEndLen"] <= self.configuration["pychic_maxFragLen"])
        ]

        if int(x.shape[0]) == 0:
            logger.fatal("All interactions have been filtered out.")

        x_shape = int(x.shape[0])

        x = x.loc[(x["distSign"] != 0) | # pylint: disable=invalid-name
                  (x["distSign"].isnull())
                 ]

        if int(x.shape[0]) < x_shape:
            logger.info("Filtered out "+
                        str(x_shape-int(x.shape[0]))+
                        " self-ligation events."
                       )

        if int(x.shape[0]) == 0:
            logger.fatal("All interactions have been filtered out.")

        ## remove baits that have no observations within the proximal range
        baits_N_sum = x.groupby("baitID").sum() # pylint: disable=invalid-name
        baits_more_min = baits_N_sum.loc[baits_N_sum["N"] >= \
                                         self.configuration["pychic_minNPerBait"]]

        x = x.loc[x["baitID"].isin(baits_more_min.index)] # pylint: disable=invalid-name


        logger.info("Filtered out "+str(baits_more_min.shape[0]-baits_N_sum.shape[0])+
                    " baits with < minNPerBait reads ("+
                    str(self.configuration["pychic_minNPerBait"])+")")

        if int(x.shape[0]) == 0:
            logger.fatal("All interactions have been filtered out.")

        # remove adjacent pairs
        if self.configuration["pychic_removeAdjacent"] is True:
            #adjacent = chinput_filt3.loc abs(chinput_filt3["baitID"]-chinput_filt3["otherEndID"])
            x = x.loc[abs(x['baitID']-x["otherEndID"]) > 1] # pylint: disable=invalid-name
            logger.info("Removed interactions with fragments adjacent to baits.")

        if int(x.shape[0]) == 0:
            logger.fatal("All interactions have been filtered out.")

        x_shape = int(x.shape[0])

        #remove baits without proximal non-bait2bait interactions
        x["isBait2bait"] = np.where(x["baitID"].isin(baitmap_id) &
                                    x["otherEndID"].isin(baitmap_id),
                                    "TRUE", "FALSE")




        NoB2BProxInter = x.loc[x["isBait2bait"] == "FALSE"] # pylint: disable=invalid-name
        NoB2BProxInter = NoB2BProxInter.loc[ # pylint: disable=invalid-name
            x["distSign"] < self.configuration["pychic_maxLBrownEst"]
        ]
        NoB2BProxInter = NoB2BProxInter[["baitID", "N"]].groupby("baitID").sum() # pylint: disable=invalid-name

        NoB2BProxInter = NoB2BProxInter[NoB2BProxInter["N"] > 0] # pylint: disable=invalid-name

        x = x[x["baitID"].isin(NoB2BProxInter.index)] # pylint: disable=invalid-name

        logger.info("Filtered out "+ str(x.shape[0]-x_shape)+
                    " baits without proximal non-Bait2bait interactions")

        if int(x.shape[0]) == 0:
            logger.fatal("All interactions have been filtered out.")

        return x

    def dataframe_merge(self, chinputs):
        """
        This function is going to execute the merge command from python to
        merge the dataframe in case there is more than onw biological replicate

        Parameter
        ---------
        Chinputs : dict
            dict with dataframes

        Return
        ------
        chunput_merged : dataframe merged
        """
        i = 0
        j = 1

        chinput_merge = pd.merge(chinputs[str(i)], chinputs[str(j)],
                                 how='outer',
                                 on=['baitID', 'otherEndID', 'otherEndLen',
                                     'distSign', 'isBait2bait'])

        if len(chinputs) == 2:
            return chinput_merge

        while j+1 < len(chinputs):
            j += 1
            chinput_merge = self.dataframe_merge_recur(
                chinput_merge,
                chinputs[str(j)])

        return chinput_merge

    @staticmethod
    def dataframe_merge_recur(chinput_merge, new_chinput):
        """
        This function is going to execute the merge command from python to
        merge the dataframe in case there is more than onw biological replicate

        Parameter
        ---------
        Chinputs : dict
            dict with dataframes

        Return
        ------
        chunput_merged : dataframe merged
        """

        chinput_merge = pd.merge(chinput_merge, new_chinput,
                                 how='outer',
                                 on=['baitID', 'otherEndID',
                                     'otherEndLen', 'distSign', 'isBait2bait'])
        return chinput_merge

    @staticmethod
    def getSampleScalingFactors(chinput_merge, pychic_maxLBrownEst, npb): # pylint: disable=invalid-name
        """
        This function calculates the scaling factors for the
        biological replicates

        Parameters
        ----------
        chinput_merge: data frame

        Returns
        -------
        s_k = list
            scaling factors for each replicate
        npb = dataframe
        """

        N_compile = re.compile(r"N.\d+") # pylint: disable=invalid-name
        N_cols = filter(N_compile.match, chinput_merge.columns) # pylint: disable=invalid-name
        ns = list(map(lambda x: x.split(".")[1], N_cols)) # pylint: disable=invalid-name
        n = len(N_cols) # pylint: disable=invalid-name

        chinput_merge = chinput_merge.ix[
            abs(chinput_merge["distSign"]) < pychic_maxLBrownEst
        ]

        column_names = ["baitID"]

        npb = pd.read_csv(npb, sep="\t", skiprows=1, header=None,
                          )

        column_names = {}
        column_names[0] = "baitID"

        for i in enumerate(npb.columns, 1):
            column_names[i[0]] = "bin"+str(i[0])

        npb.rename(columns=column_names, inplace=True)

        N_p_bait = pd.DataFrame() # pylint: disable=invalid-name

        npb["sum"] = npb.iloc[:, 1:].sum(axis=1)

        N_p_bait = npb[["baitID", "sum"]] # pylint: disable=invalid-name

        npb.drop("sum", 1, inplace=True) # pylint: disable=invalid-name

        chinput_merge = pd.merge(chinput_merge, N_p_bait, how="left",
                                 on="baitID")

        N_p_bait = chinput_merge[["baitID", "sum"]+N_cols] # pylint: disable=invalid-name

        N_p_bait = N_p_bait.groupby(["baitID", "sum"], # pylint: disable=invalid-name
                                    as_index=False)[N_cols].sum()

        #Now divide the total number of contacts per Baits by the total OE per baits
        for i in range(1, n+1): # pylint: disable=invalid-name
            N_p_bait["s_"+str(i)+"j"] = N_p_bait["N."+str(i)]/N_p_bait["sum"]

        s_Nj = ["s_"+i.split(".")[1]+"j" for i in N_cols] # pylint: disable=invalid-name

        N_p_bait["geo_mean"] = stats.gmean(N_p_bait[s_Nj], axis=1)

        s_k = {}

        for i in enumerate(ns, 1):
            array_m = np.array(N_p_bait["s_"+str(i[0])+"j"]/N_p_bait["geo_mean"])
            #Division 0/0 generated inf element, which makes the median sligthly different
            #from R, thats why im extracting the finite numbers
            condition = np.isfinite(array_m) is True
            array_m = np.extract(condition, array_m)
            s_k["N."+str(i[0])] = np.median(array_m)

        return s_k, npb

    def merge_chinputs(self, chinputs, npb):
        """
        This function will work in case there is more than one chinput file.
        It takes the filtered chinput files and merge them

        Parameters
        ----------
        chinputs: list

        Returns
        -------
        chinputs_merged: data frame
        """
        for i in enumerate(chinputs):
            column_names = chinputs[str(i[0])].columns
            #check if the column N is already N.<k>. k being the i
            if "N" not in column_names:
                N_compile = re.compile(r"N.\d+") # pylint: disable=invalid-name
                N_name = filter(N_compile.match, column_names) # pylint: disable=invalid-name
                if not N_name:
                    logger.fatal("Did not found 'N' coulmn in chinput file")
            else:
                chinputs[str(i[0])] = chinputs[str(i[0])].rename(columns={"N" : "N."+str(i[0]+1)})

            # To merge the two dataframes, first convert NaNs values in "NA".
            # Cause we dont want to mix the NaN values the the dataframes already
            # have, with the NaN generated in the merge execution. Then merge the
            # two dataframes, then convert NaNs into 0, and "NA" into NaN.

            chinputs[str(i[0])] = chinputs[str(i[0])].replace(np.nan, "NA")

        chinput_merge = self.dataframe_merge(chinputs)

        #replace NaN with 0, and "NA" with NaN
        chinput_merge = chinput_merge.replace(np.nan, 0)
        chinput_merge = chinput_merge.replace("NA", np.nan)

        logger.info("computing merged scores..")

        s_k, npb = self.getSampleScalingFactors(chinput_merge,
                                                self.configuration["pychic_maxLBrownEst"],
                                                npb)

        #New column with notmalize merged counts
        #N.1*s_ks["N.1"]+N.2*s_ks["N.2"]+N.3*s_ks["N.3"])/sum(s_ks))
        N_compile = re.compile(r"N.\d+") # pylint: disable=invalid-name
        N_name = filter(N_compile.match, chinput_merge.columns) # pylint: disable=invalid-name

        #create a N column with the weighted mean
        # N = round(N.1*s_ks["N.1"]+N.2*s_ks["N.2"]+N.3*s_ks["N.3"])/sum(s_ks))
        for i in enumerate(s_k, 1):
            chinput_merge["N."+str(i[0])] = chinput_merge["N."+str(i[0])]*s_k["N."+str(i[0])]


        chinput_merge["N"] = chinput_merge[N_name].sum(axis=1)/sum([i for i in s_k.values()])
        chinput_merge["N"] = chinput_merge["N"].round()
        chinput_merge["N"] = np.where(chinput_merge.N == 0, 1, chinput_merge.N)

        return chinput_merge, npb

    def normaliseFragmentSets(self, x, viewpoint, idcol, Ncol, binsize, # pylint: disable=invalid-name
                              npb=False, adjBait2bait=True,
                              shrink=False, refExcludeSuffix=None,
                              plot=True, outfile=None,
                              debug=False):
        """
        The normalisation engine used for normaliseBaits and normaliseOtherEnds
        "Viewpoint" will be used in the comments for either baits or sets of other ends,
        depending on the direction of normalisation.

        npb is a data table containing the number of reads per viewpoint per bin,
        for viewpoint=="bait" it's just the table from the nperbinfile (read via .readNPBfile),
        where npb's idcol should match x's idcol.
        for viewpoit=="otherEnd" the table from the nbaitsperbin file needs preprocessing
        to sum over pools of other ends, and this is done in normaliseOtherEnds,
        with the nbp column added to x before submitting to .normaliseFragmentSets.

        s is the current chicagoData object's settings list

        Parameters
        ----------
        x : dataframe
        viewpoint : str
            bait or OtherEnd
        npb : str
            path to npb file
        idcol: str
            baitID otherEndID

        Returns
        -------
        x_
        """
        if viewpoint == "bait":
            scol = "s_j"

            x["distbin"] = pd.cut(x["distSign"].abs(),
                np.arange(0, self.configuration["pychic_maxLBrownEst"]+1, binsize))

            x["distbin"] =  x['distbin'].apply(lambda x: x.left)

            x["distbin"] = x["distbin"].astype(float)

            xAll = x # pylint: disable=invalid-name

            if adjBait2bait:
                """
                BE CAREFULL IF THIS FUNCTION IS CALLED OUTSIDE NORMALISEDBAITS
                if "isBait2bait" not in x.columns:
                    x[, isBait2bait := FALSE]
                    x[wb2b(otherEndID), isBait2bait:= TRUE]
                """
                x = x.loc[x["isBait2bait"] == "FALSE"]

            # x is the data table used to compute the scaling factors
            x = x.loc[x["distbin"].notna()]

            #x["bincol"] = ["bin"+str(int(i/binsize)+1) for i in x["distbin"]]

            x["bincol"] = [int(i/binsize)+1 for i in x["distbin"]]

            npb_dic = {}
            with open(npb, "r") as npb_file:
                for line in npb_file:
                    line_hdl = line.rstrip().split("\t")
                    try:
                        npb_dic[int(line_hdl[0])] = line_hdl[1:]
                    except ValueError:
                        continue

            ntot = []

            for r in zip(x["baitID"], x["bincol"]):
                ntot.append(npb_dic[r[0]][r[1]-1])

            x["ntot"] = [int(i) for i in ntot]

        else:
            scol = "s_i"

        logger.info("Computing binwise means..")

        x = x[[idcol, Ncol, "distbin", "ntot"]]

        #sbbm is the input matrix for normalisation that contains the mean number of reads per bin for each bait
        #This makes the dataframe sum all the interactions "N" ma
        sbbm = x.groupby(
            ["ntot", "distbin", idcol],
            sort=False,
            as_index=False).sum()

        sbbm["bbm"] = sbbm[Ncol]/sbbm["ntot"]

        sbbm.drop([Ncol, "ntot"], axis=1, inplace=True)

        sbbm.sort_values(by=["distbin"], inplace=True)

        if viewpoint == "bait" or refExcludeSuffix == False:
            #IMPROVE THIS APPLYING THE FUNCTION TO BLOCK AND CREATING A COLUMN AT THE SAME TIME
            geomean = sbbm.groupby(
                ["distbin"],
                sort=False,
                as_index=False).bbm.apply(gmean)

            bins_elements = Counter(sbbm["distbin"])

            geo_mean = []
            m = 0
            for i in range(0, int(max(bins_elements))+1, binsize):
                geo_mean = geo_mean + [geomean[m]]*bins_elements[i]
                m +=1

            sbbm["geo_mean"] = geo_mean

        else:

            distbin_sorted = sorted(sbbm["distbin"].unique())

            sbbm_NB2B = sbbm[sbbm["tlb"].astype("str").str.contains("B2B")==False]

            geomean = sbbm_NB2B.groupby(
                ["distbin"],
                sort=False,
                as_index=False)["bbm"].apply(gmean)

            geomean = pd.DataFrame(geomean, columns=["geo_mean"])

            geomean["distbin"] = distbin_sorted

            sbbm = pd.merge(sbbm, geomean, how="left", on="distbin")

        #DEseq-style normalisation
        if not shrink or viewpoint == "otherEnd":
            sbbm["s_iv"] = sbbm["bbm"]/sbbm["geo_mean"]
            s_v = sbbm.groupby(idcol, as_index=False).s_iv.median()
        else:
            logger.info("computing shrunken means...")
            ##
            #
            #
            #
            #
            #.....

        s_v.rename(columns={"s_iv": scol}, inplace=True)

        if len(s_v[s_v[scol].isnull()]) > 0:
            logger.info("The following viewpoints couldn't be robustly "
                        "normalised (too sparse?) and will be removed:")
            logger.info(s_v[s_v[scol].isnull()])
            s_v = s_v[s_v[scol].notnull()]

        if viewpoint == "otherEnd":
            xAll = x
            xAll = pd.merge(xAll, s_v, how="left", on="tlb")

        if viewpoint == "bait":
            xAll = pd.merge(s_v, xAll, how="left",
                            on = "baitID")
            gm = pd.DataFrame({"refBinMean" : geomean,
                               "distbin" : pd.unique(sbbm["distbin"])})
            xAll = pd.merge(xAll, gm, how="left", on="distbin")


        return xAll


    def normaliseBaits(self, x, npb):
        """
        This function normalize the baits, calulating the s_j parameter part of the
        mean of the negative binomial model used to model the background error

        Parameter
        ---------
        x: str
            path to the chinput file already filtered

        Return
        ------
        chinput_
        """
        binsize = self.configuration["pychic_binsize"]

        x = self.normaliseFragmentSets(x, "bait", "baitID",
                                       "N", binsize, npb)

        x["NNb"] = (x["N"]/x["s_j"]).round()
        x["NNb"] = np.where(x["NNb"] == 0, 1, x["NNb"])

        return x


    def addTLB(self, chinput_j, adjBait2bait=True):
        """
        ##Assigns each fragment a "tlb" - a range containing the number of trans
        ##1) The bins are constructed based on reduced data - (outliers trimmed, )
        ##2) The bin endpoints are readjusted such that no fragments fall outside.
        ##3) These bins are then applied to the entire dataset.

        Parameters
        ----------
        chinput_j: dataFrame
        """
        filterTopPercent = self.configuration["pychic_tlb_filterTopPercent"]
        minProxOEPerBin = self.configuration["pychic_tlb_minProxOEPerBin"]
        minProxB2BPerBin = self.configuration["pychic_tlb_minProxB2BPerBin"]

        logger.info("preprocessing input....")

        # Checking whether in the input, we had distances at trans-interactions labeled as NA
        # (as opposed to a dummy maximum distance)
        transNA = False

        if chinput_j["distSign"].isnull().values.any():
            transNA = True
            transD = max(chinput_j["distSign"])+self.configuration["pychic_binsize"]
            chinput_j.distSign.fillna(transD, inplace=True)
            #chinput_j["distSign"] = np.where(chinput_j["distSign"].isna(), transD)


        else:
            transD = max(chinput_j["distSign"])
            logger.info("Warning: No NAs found in input. "
                        "Assuming the max distance of "+str(transD)+
                        " is a dummy for  trans-counts.")

        logger.info("Computing trans-counts")

        if adjBait2bait:
            if not "isBait2bait" in chinput_j.columns:
                #chinput_j = chinput_j.ix[chinput_j["is"]]
                #wb2b
                pass

            chinput_j["transLength"] = np.where(chinput_j['distSign'] == transD, 1, 0)

            transLen = chinput_j[["otherEndID",
                                  "transLength",
                                 ]].groupby(
                                     "otherEndID", as_index=False).sum()

            transLen2 = pd.DataFrame({"otherEndID" : chinput_j["otherEndID"], \
                                     "isBait2bait" : chinput_j["isBait2bait"], \
                                     "distSign" : abs(chinput_j["distSign"])})


            transLen2.sort_values(["otherEndID","distSign"], inplace=True)

            transLen2 = transLen2[["distSign" , "isBait2bait", "otherEndID"]]

            distSign = []
            isBait2bait = []
            otherEndID = []
            counter= 0
            for row in transLen2.itertuples():
                counter += 1
                if counter == 1 :
                    distSign.append(row[1])
                    isBait2bait.append(row[2])
                    otherEndID.append(row[3])
                    OE = row[3]
                else:
                    if row[3] == OE:
                        continue
                    else:
                        distSign.append(row[1])
                        isBait2bait.append(row[2])
                        otherEndID.append(row[3])
                        OE = row[3]

            transLen2 = pd.DataFrame({"otherEndID" : otherEndID,
                                      "isBait2bait" : isBait2bait,
                                      "distSign" : distSign
                                     })

            transLen = pd.merge(transLen, transLen2, on="otherEndID", how="right")


        else:
            pass
            #something similar to the above function

        #remember to rename the transD column for te else statement
        transLend0 = transLen.copy()

        transLen = transLen[transLen["transLength"] <= np.percentile(
            transLen["transLength"], 100-self.configuration["pychic_tlb_filterTopPercent"])]


        if adjBait2bait:
            transLend0 = transLen.copy()
            transLen = transLend0.loc[transLend0["isBait2bait"] == "FALSE"]
            transLenB2B = transLend0[transLend0["isBait2bait"] == "TRUE"]

        logger.info("Binning..")

        distSign = transLen[transLen["distSign"].abs() <= self.configuration["pychic_maxLBrownEst"]]["transLength"]

        cuts = self.cut2(distSign, int(self.configuration["pychic_tlb_minProxOEPerBin"]))

        #with depleated datasets curs == 1
        if len(cuts) == 1:
            tlbClasses = cuts*transLen.shape[0]
        else:
            min_transLength = min(transLen["transLength"])
            max_transLength = max(transLen["transLength"])
            if min(cuts) > min_transLength:
                cuts[0] = min_transLength
            if max(cuts) < max_transLength:
                cuts[-1] = max_transLength
            tlbClasses = pd.cut(transLen["transLength"],
                                cuts,
                                include_lowest=True
                                )

        transLen.loc[:, "tlb"] = tlbClasses


        if adjBait2bait:

            # cutB2B = cut2(transLenB2B[abs(distSign)<= s$maxLBrownEst]$transLength,
            #        m=minProxB2BPerBin, onlycuts=TRUE)
            #cutB2B = [0, 1, 2, 3, 4, 5, 6, 7, 8]

            transLen_length = transLenB2B[transLenB2B["distSign"].abs() <= \
                self.configuration["pychic_maxLBrownEst"]]["transLength"]

            cutB2B = self.cut2(transLen_length, self.configuration["pychic_tlb_minProxB2BPerBin"] )

            if len(cutB2B) == 1:
                tlbClassesB2B = cutB2B*transLenB2B.shape[0]
            else:
                min_transLength = min(transLenB2B["transLength"])
                max_transLength = max(transLenB2B["transLength"])
                if min(cutB2B) > min_transLength:
                    cutB2B[0] = min_transLength
                if max(cutB2B) < max_transLength:
                    cutB2B[-1] = max_transLength

                tlbClassesB2B = pd.cut(transLenB2B["transLength"],
                                       cutB2B,
                                       include_lowest=True)

            transLenB2B.loc[:, "tlb"] = tlbClassesB2B

            transLenB2B.loc[:, "tlb"] = transLenB2B["tlb"].astype("str") + "B2B"

            transLen = transLen.append(transLenB2B)

        transLen.drop(["transLength", "isBait2bait", "distSign"],
                      axis=1,
                      inplace=True
                     )

        chinput_j.drop(["transLength"], axis=1, inplace=True)

        # Discard TLB column if present in x
        #if "tlb" in chinput_j delete it
        if "tbl" in chinput_j.columns:
            chinput_j.drop(["tbl"], axis=1, inplace=True)


        chinput_j = pd.merge(chinput_j, transLen, how="left", on="otherEndID")


        """
         428139 0.9398174     419466 1        3195 -31823781       FALSE         NA
        2  427751 1.2124222     419466 1        3195 -30521967       FALSE         NA
        3  425079 2.0696748     419466 1        3195 -21272923       FALSE         NA
        4  424854 1.2550334     419466 2        3195 -20594316       FALSE         NA
        5  423617 0.6890949     419466 1        3195 -16816825       FALSE         NA
        6  423560 1.2931352     419466 1        3195 -16591004       FALSE         NA
        7  423450 1.3431198     419466 1        3195 -16216369       FALSE         NA
        8  421902 2.6059243     419466 1        3195 -11008896       FALSE         NA
        9  403965 1.3534728     419466 1        3195  62753699       FALSE         NA
        10 404491 1.6695769     419466 1        3195  62753699       FALSE         NA
        11 408471 0.5752988     419466 1        3195  62753699       FALSE         NA
        12 409144 0.8654104     419466 1        3195  62753699       FALSE         NA
        13 410108 0.6957737     419466 1        3195  62753699       FALSE         NA
        14 411133 0.5165262     419466 1        3195  62753699       FALSE         NA
        15 411164 0.9946355     419466 1        3195  62753699       FALSE         NA
        16 411831 1.5086577     419466 1        3195  62753699       FALSE         NA
        17 412613 2.0638962     419466 1        3195  62753699       FALSE         NA
        18 418022 0.7670651     419466 1        3195  62753699       FALSE         NA
        19 418110 0.7582941     419466 1        3195  62753699       FALSE         NA
        """
        chinput_j =  chinput_j[chinput_j["otherEndID"] != 419466]

        if transNA:
            chinput_j["distSign"] = np.where(chinput_j["distSign"] == max(chinput_j["distSign"]), \
                                             np.nan, \
                                             chinput_j["distSign"])
        return chinput_j


    def normaliseOtherEnds(self, chinput_j, nbpb, Ncol="NNb", normNcol="NNboe"):
        """
        This function is used to calculate the normalization parameter for
        the OE

        Parameters
        ----------
        chinput_j: dataframe
            dataframe after normalize the Baits
        Ncol: str
            name of the Ncol colum
        normNcol: str
            name of the normalized column

        Returns
        -------
        DataFrame
        """

        chinput_j = self.addTLB(chinput_j)

        x = chinput_j.ix[(chinput_j["distSign"].abs() <= self.configuration["pychic_maxLBrownEst"]) &
                         (chinput_j["distSign"].notnull())
                        ]

        logger.info("Computing total bait counts...")

        nbpb = self.prepare_design(nbpb, "otherEndID")

        nbpb = pd.merge(x, nbpb, how="left", on="otherEndID")

        # Compute the sums of observed baits per bin for pools of other ends -
        # NB: we need to some only once for each other end, but each other end is present
        # more than once for each tlb

        #Eliminate the first occurence from OE and distbin
        nbpb.sort_values(["otherEndID", "distbin"], ascending =[True, True], inplace =True)

        nbpb.drop_duplicates(subset = ['otherEndID', 'distbin'], keep="first", inplace=True)

        nbpb["BinN"] = nbpb["distbin"]/self.configuration["pychic_binsize"]+1

        nbpb["BinN"] = "bin"+ nbpb["BinN"].astype(int).astype(str)

        nbpb.reset_index(inplace=True)


        ntot = []
        for binN in nbpb["BinN"].iteritems():
            ntot.append(nbpb.loc[binN[0],binN[1]])

        nbpb["ntot"] = ntot

        nbpbSum = nbpb[["tlb","distbin","ntot"]]

        nbpb.drop(["ntot"], axis=1, inplace=True)

        nbpbSum = nbpbSum.groupby(["tlb","distbin"], as_index=False)["ntot"].sum()

        x = pd.merge(x, nbpbSum, how="left", on=["tlb","distbin"])

        logger.info("Computing scaling factors")

        x =  self.normaliseFragmentSets(x, "otherEnd", "tlb", "NNb", self.configuration["pychic_binsize"],
                                        refExcludeSuffix="B2B")

        x.drop_duplicates(subset=["s_i"], keep="first", inplace=True)

        x.drop(["NNb", "distbin", "ntot"], axis=1, inplace=True)

        #PLOT
        logger.info("Computing normalised counts...")

        if "s_i" in chinput_j.columns:
            chinput_j.drop(["s_i"], axis=1, inplace=True)


        chinput_j = pd.merge(chinput_j, x, how="left", on="tlb")

        #If we can estimate s_i robustly, assume it to be one
        chinput_j["tlb"].fillna(1, inplace=True)


        chinput_j["NNboe"] = chinput_j[Ncol]/chinput_j["s_i"]
        chinput_j["NNboe"]  = chinput_j["NNboe"].round(decimals=0)
        chinput_j["NNboe"] = np.where(chinput_j["NNboe"] > 1, chinput_j["NNboe"], 1 )

        return chinput_j


    def prepare_design(self, design, first_col):
        """
        This function is runned in case that there is just one chinput.
        Becuase the npb file is modified during the merge_input function

        Parameters
        ----------
        npb: str
            path to the npb file
        Returns
        -------
        DataFrame
        """
        design = pd.read_csv(design, sep="\t", skiprows=1, header=None,
                         )

        column_names = {}
        column_names[0] = first_col

        for i in enumerate(design.columns,1):
            column_names[i[0]] = "bin"+str(i[0])

        design.rename(columns=column_names, inplace=True)

        return design

    def estimateTechnicalNoise(self, chinput_ji, rmap, baitmap):
        """
        This function estimate the technical noise of the experiments

        Parameters
        ----------
        chinput_ji : DataFrame
            output from normaliseOtherEnds function

        Returns
        -------
        chinput_jiw

        """
        logger.info("Estimating technical noise based on trans-counts...")

        minBaitsPerBin = self.configuration["pychic_techNoise_minBaitsPerBin"]
        adjBait2bait = self.configuration["pychic_adjBait2bait"]

        logger.info("Binning baits based on observed trans-counts...")


        trans = chinput_ji[chinput_ji["distSign"].isnull()]

        transBaits = trans['baitID'].value_counts()

        transBaitLen = pd.DataFrame(transBaits)

        transBaitLen = transBaitLen.rename(columns = {"baitID": "transBaitLen"})

        transBaitLen["baitID"] = transBaitLen.index.tolist()

        #levels = [6, 22, 32, 46, 131]
        levels = self.cut2(transBaitLen["transBaitLen"],
                           self.configuration["pychic_techNoise_minBaitsPerBin"])

        transBaitLen["tblb"] = pd.cut(transBaitLen["transBaitLen"], levels, right=False)
        transBaitLen["tblb"] = transBaitLen["tblb"].astype(str)
        transBaitLen["tblb"] = np.where(transBaitLen["tblb"] == "nan",
                                        #"[46, 131)"
                                        "["+str(levels[-2])+", "+str(levels[-1])+")",
                                        transBaitLen["tblb"])
        chinput_ji = pd.merge(chinput_ji, transBaitLen, how="left", on="baitID")

        logger.info("Defining interaction pools and gathering "+
                    "the observed numbers of trans-counts per pool...")

        #Getting the observed numbers of trans-counts
        #Mistery of Ntrans

        logger.info("Computing the total number of possible interactions per pool...")
        logger.info("Preparing the data...")

        baitmap = pd.read_csv(baitmap, sep="\t", header=None)
        rmap = pd.read_csv(rmap, sep="\t", header=None)

        baitmap = baitmap.rename(columns = {0:"baitchr", 3:"baitID"})
        baitmap.drop([1, 2, 4], axis=1, inplace=True)
        rmap = rmap.rename(columns = {0:"otherEndchr", 3:"otherEndID"})
        rmap.drop([1, 2], axis=1, inplace=True)

        chinput_ji = pd.merge(chinput_ji, rmap, how="left", on="otherEndID")
        chinput_ji = pd.merge(chinput_ji, baitmap, how="left", on="baitID")

        logger.info("\nProcessing fragment pools")

        chinput_ji["tlb"] = chinput_ji["tlb"].astype(str)

        res_chinput = chinput_ji[chinput_ji["distSign"].isnull()]
        res = res_chinput.groupby(["tlb", "tblb"], as_index=False)["N"].sum()

        tlb_tblb = chinput_ji.drop_duplicates(["tlb", "tblb"])
        tlb_tblb = tlb_tblb[["tlb","tblb"]]


        numPairsdf = pd.DataFrame(columns = {"tlb", "tblb", "numPairs"})

        for i in tlb_tblb.index:
            tlb = tlb_tblb.at[i, "tlb"]
            tblb = tlb_tblb.at[i, "tblb"]
            temp_chinput = chinput_ji[(chinput_ji["tlb"] == tlb) &
                                      (chinput_ji["tblb"] == tblb)]

            baits = temp_chinput["baitID"].unique()
            oes = temp_chinput["otherEndID"].unique()

            baits_set = set(baits)
            oes_set = set(oes)
            bChr = temp_chinput.drop_duplicates("baitID")
            bChr = bChr["baitchr"]
            oeChr = temp_chinput.drop_duplicates("otherEndID")
            oeChr = oeChr["otherEndchr"]

            numPairs = 0
            for chromo in bChr.unique():
                pairs = (len([x for x in bChr if x== chromo])*len([y for y in oeChr if y != chromo]))-len(baits_set.intersection(oes_set))
                numPairs += pairs

            numPairsdf = numPairsdf.append({
                "tlb" : tlb,
                "tblb" : tblb,
                "numPairs" : numPairs
                }, ignore_index=True)

        res = pd.merge(numPairsdf, res,  how="left", on=["tlb", "tblb"])

        res = res.rename(columns={"N" : "nTrans"})

        res["Tmean"] = res["nTrans"]/res["numPairs"]

        res.drop(["numPairs",
                  "nTrans"], axis=1, inplace=True)

        #PLOT
        logger.info("Post-processing the results...")

        chinput_ji = pd.merge(chinput_ji, res, how="left", on=["tlb", "tblb"])

        chinput_ji.drop(["transBaitLen",
                         "baitchr",
                         "otherEndchr"],
                         axis=1,
                         inplace=True)

        return chinput_ji

    def estimateDistFun(self, chinput_jiw):
        """
        Estimate the distancde function

        Parameters
        ----------
        chinput_jiw: DataFrame
            output from estimateTechnicalNoise
        """
        # Take the "refBinMean" column of the data x as f(d_b)
        # then interpolate & extrapolate to get f(d).
        # TODO output extra diagnostic information?

        # Get f(d_b)
        chinput_jiw = chinput_jiw[chinput_jiw["refBinMean"].notnull()]
        f_d = chinput_jiw.drop_duplicates(["distbin", "refBinMean"])
        f_d = f_d[["distbin", "refBinMean"]]
        f_d.sort_values(by=["refBinMean"], inplace=True, ascending=False)

        f_d["midpoint"] = np.arange(self.configuration["pychic_binsize"]/2,
                             self.configuration["pychic_binsize"]*f_d.shape[0],
                             self.configuration["pychic_binsize"])

        obs_min = np.log(min(f_d["midpoint"]))
        obs_max = np.log(max(f_d["midpoint"]))

        distFunParams = {}
        ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
        f_d_cubic = np.polyfit(np.log(f_d["midpoint"]), np.log(f_d["refBinMean"]), 3)
        fit = [f_d_cubic[3], f_d_cubic[2], f_d_cubic[1], f_d_cubic[0] ]
        distFunParams["cubicFit"] = fit

        distFunParams["obs_min"] = obs_min
        distFunParams["obs_max"] = obs_max

        beta1 = fit[1] + 2*fit[2]*obs_min + 3*fit[3]*(obs_min**2)
        beta2 = fit[1] + 2*fit[2]*obs_max + 3*fit[3]*(obs_max**2)

        alpha1 = fit[0] +(fit[1] - beta1)*obs_min + fit[2]*obs_min**2 +fit[3]*obs_min**3
        alpha2 = fit[0] +(fit[1] - beta2)*obs_max + fit[2]*obs_max**2 +fit[3]*obs_max**3

        distFunParams["head_coef"] = [alpha1, beta1]
        distFunParams["tail_coef"] = [alpha2, beta2]

        #PLOT

        return distFunParams


    def readProxOEfile(self, poe):

        """
        Reads a pre-computed text file that denotes which other ends are in the proximal
        range relative to each bait, and gives that distance.
        Note that fragments that are too small/too large have already been removed.
        s is the current chicagoData object's settings list

        Parameters
        ----------
        poe : str
            poe file generated my makeDesingFiles.py

        Returns
        -------
        poe : DataFrame
        """
        with open(poe, "r") as file_poe:
            params = file_poe.readline()

        params = [param for param in params.split("\t")]

        minsize = "".join([param for param in params if param.startswith("minFragLen")])
        if int(minsize.split("=")[1]) != self.configuration["pychic_minFragLen"]:
            logger.info("The minFragLen specified in the ProxOE file header is not "+
                        "equal to minFragLen defined in experiment settings. Amend either "+
                        "parameter setting (and if needed, generate a new ProxOE file) before"+
                        "running the analysis\n")

        maxsize = "".join([param for param in params if param.startswith("maxFragLen")])
        if int(maxsize.split("=")[1]) != self.configuration["pychic_maxFragLen"]:
            logger.info("The maxFragLen specified in the ProxOE file header is not"+
                        " equal to maxFragLen defined in experiment settings."+
                        " Amend either parameter setting (and if needed, generate"+
                        " a new ProxOE file) before running the analysis\n")

        maxl = "".join([param for param in params if param.startswith("maxLBrownEst")])
        if int(maxl.split("=")[1]) != self.configuration["pychic_maxLBrownEst"]:
            logger.info("The maxLBrownEst specified in the ProxOE file header is "+
                        "not equal to maxLBrownEst defined in experiment settings."+
                        "Amend either parameter setting (and if needed, generate a "+
                        "new ProxOE file) before running the analysis\n")

        binsz = "".join([param for param in params if param.startswith("binsize")])
        if int(binsz.split("=")[1]) != self.configuration["pychic_binsize"]:
            logger.info("The binsize specified in the ProxOE file header is"+
                       "not equal to binsize defined in experiment settigs."+
                       " Amend either parameter setting (and if needed, "+
                       "generate a new ProxOE file) before running the analysis\n")

        removeb2b = "".join([param for param in params if param.startswith("removeb2b")])
        if removeb2b.split("=")[1] != "True":
            logger.info("The ProxOE file must be generated with removeb2b==True."+
                        "Please generate a new file.\n")

        removeAdjacent = "".join([param for param in params if param.startswith("removeb2b")])
        if removeAdjacent.split("=")[1] == True and self.configuration["pychic_removeAdjacent"] is False:
            logger.info("The removeAdjacent parameter settings used for generating ProxOE " +
                        "file (according to its header) and defined in experiment settings "+
                        "do not match. Amend either setting (and if needed, generate a new "+
                        "ProxOE file) before running the analysis\n")


        rmapfile = "".join([param for param in params if param.startswith("rmapfile")])
        if os.path.split(rmapfile.split("=")[1])[1] != os.path.split(input_files["RMAP"])[1]:
            logger.info("The .rmap files used for generating the ProxOE file "+
                        "(according to the ProxOE header) and the one defined "+
                        "in experiment settings do not match. Amend either setting "+
                        "(and if needed, generate a new ProxOE file) before running the analysis\n")

        poe = pd.read_csv(poe, sep="\t", skiprows=1, header=None, names=["baitID", "otherEndID", "dist"])

        return(poe)

    def distFun(self, d, distFunParams):
        """
        Select the distance function parameters

        Parameters
        ----------
        disFunParams: dic
        d: series
            x["distSign"] column

        Returns
        -------
        out : list
        """
        obs_max = distFunParams["obs_max"]
        obs_min = distFunParams["obs_min"]
        head_coef = distFunParams["head_coef"]
        tail_coef = distFunParams["tail_coef"]
        fit = distFunParams["cubicFit"]

        ##Put everything together to get the final function
        d = np.log(d)

        out = []

        for dist in d:
            if dist > obs_max:
                out.append(tail_coef[0] + dist*tail_coef[1])
            elif dist < obs_min:
                out.append(head_coef[0] + dist*head_coef[1])
            else:
                out.append(fit[0] + fit[1]*dist + fit[2]*(dist**2) + fit[3]*(dist**3))


        return np.exp(out)


    def estimateBMean(self, x, distFunParams):
        """
        Estimate the Bronian Mean

        Parameters
        ----------
        x: DataFrame

        Returns
        -------
        x: Dataframe
        """

        if "s_i" in x.columns:
            out = self.distFun(x["distSign"].abs(), distFunParams)
            x["Bmean"] = x["s_j"]*x["s_i"]*out
        else:
            logger.info("s_i factors NOT found in .estimateBMean - variance will increase, estimating means anyway...")
            out = self.distFun(x["distSign"], distFunParams)
            x["Bmean"] = x["s_j"]*out

        x["Bmean"] = x["Bmean"].fillna(0)

        return x

    def estimateDispersion(self, chinput_jiw, proxOE, distFunParams):
        """
        Estimate the dispersion

        Parameters
        ----------
        chinput_jiw: DataFrame
        poe: DataFrame

        Returns
        -------
        params
        """

        if "s_i" in chinput_jiw.columns:
            siPresent = True
        else:
            siPresent = False

        adjBait2bait = self.configuration["pychic_adjBait2bait"]
        samples = self.configuration["pychic_brownianNoise_samples"]
        subset = self.configuration["pychic_brownianNoise_subset"]
        subset = False
        maxLBrownEst = self.configuration["pychic_maxLBrownEst"]

        ##Pre-filtering: get subset of data, store as x
        ##---------------------------------------------

        if subset != False:
            if len(chinput_jiw["baitID"].unique()) > subset:
                sel_sub = sorted(sample(chinput_jiw["baitID"].unique(), subset))
                x = chinput_jiw[chinput_jiw["baitID"].isin(sel_sub)]
            else:
                x = chinput_jiw
                subset = None
        else:
            x = chinput_jiw

        ##consider proximal region only...
        x = x[(x["distSign"].abs() < self.configuration["pychic_maxLBrownEst"]) & \
              (x["distSign"].notnull())]

        #remove bait2bait...
        if adjBait2bait == True:
            if "isBait2bait" not in x.columns:
                x = x[x["isBait2bait"] == "FALSE"]
                #CARFEULL IF THIS IS USED IN A DIFFERENT OCASION
            else:
                x = x[x["isBait2bait"] == "FALSE"]


        ##1) Reinstate zeros:
        ##----------------
        ##1A) Choose some (uncensored) baits. Pick relevant proxOE rows. Note: censored fragments,
        ##   censored bait2bait pairs (etc...) already taken care of in pre-computation of ProxOE.

        if subset != False:
            ## if we chose a subset of baits, restrict to that (none of these should be censored)
            sel_baits = sel_sub
        else:
            sel_baits = x["baitID"].unique()

        proxOE = proxOE[proxOE["baitID"].isin(x["baitID"])]

        x.sort_values(inplace=True, by=["baitID", "otherEndID"])

        sjLookup = x[["baitID", "s_j"]]
        sjLookup = sjLookup.drop_duplicates(["baitID"])

        if siPresent:
            siLookup = x[["otherEndID", "s_i"]]
            siLookup = siLookup.drop_duplicates(["otherEndID"])

        if siPresent:
            x = x[["baitID","otherEndID", "s_i", "s_j", "N", "distSign"]]
            x = pd.merge(proxOE, x, how="left", on=["baitID", "otherEndID"])
        else:
            x = x[["baitID","otherEndID", "s_j", "N", "distSign"]]
            x = pd.merge(proxOE, x, how="left", on=["baitID", "otherEndID"])

        ##Merging like this means that we are missing N, s_i, s_j information for most of the rows. So:
        ##1C) Repopulate table with information...

        # TODO: Recast following, avoid lookup tables, instead subset and assign by reference
        ## - 0s in Ncol

        x["N"] = x["N"].fillna(0)

        #siLookup.rename(inplace=True, columns={"s_1": "s_i2"})
        x.drop(["s_j", "s_i"], axis=1, inplace=True)
        x = pd.merge(x, sjLookup, how="left", on=["baitID"])

        if siPresent:
            x = pd.merge(x, siLookup, how="left", on=["otherEndID"])
            if x["s_i"].isnull().values.any():
                x["s_i"] = x["s_i"].fillna(1)

        ## - distances
        ##Sanity check - the distances should agree (modulo rounding)
        distances = x["dist"] - x["distSign"].abs()
        distances.dropna(inplace=True)
        distances = distances.abs()

        if (distances.values > 1).any():
            logger.info("estimateBrownianComponent: Distances in precomputed ProxOE file did not match distances supplied.")

        x["distSign"] = x["dist"]

        x.drop(inplace=True, columns=["dist"], axis=1)
        ##2) Calculate Bmeans
        ##----------------

        #to debug
        x.sort_values(["baitID", "otherEndID"], inplace=True)

        x = self.estimateBMean(x, distFunParams)
        ##3)Fit model
        ##---------
        x["Bmeanlog"] = np.log(x["Bmean"])

        x["zeros"] = np.zeros(len(x["Bmean"]))

        logger.info("Sampling the dispersion...")
        """
        LOOK AT THE TYPES OF THE PANDAS COLUMNS AND SEE IF THAT AFFECTS
        from rpy2.robjects.packages import importr
        from rpy2.robjects.numpy2ri import numpy2ri
        from rpy2 import robjects

        robjects.numpy2ri.activate()

        MASS = importr('MASS')
        stats = importr('stats')

        from rpy2.robjects import r
        numpyN = np.array(x["N"])
        numpyBmenalog = np.array(x["Bmeanlog"])

        r_y = numpy2ri(numpyN)
        r_x = numpy2ri(numpyBmenalog)

        r.assign("y", r_y)
        r.assign("x", r_x)

        r("x <- as.matrix(x)")
        r("y <- as.matrix(y)")
        r("res <- glm.nb(formula=y~x + 0)")
        model_theta = r("res$theta")

        print(model_theta)
        """
        model_theta = 2.5563913

        return model_theta

    def estimateBrownianComponent(self, chinput_jiw, distFunParams):
        """
        1) Reinstate zeros
        2) Add a "Bmean" column to x, giving expected Brownian component.
        3) Calculate dispersion by regressing against "Bmean", added to x as "dispersion" attribute
        subset: Since we don't need the entire data set, can just calculate
        based on a random subset of baits.
        !!NB!! Use set.seed to force subset analysis to be reproducible

        Parameters
        ----------
        chinput_jiw: DataFrame
            output of estimateTechnicalNoise

        Returns
        -------
        chinput_jiwb
        """

        adjBait2bait = self.configuration["pychic_adjBait2bait"]
        samples = self.configuration["pychic_brownianNoise_samples"]
        subset = self.configuration["pychic_brownianNoise_subset"]
        seed = self.configuration["pychic_brownianNoise_seed"]
        maxLBrownEst = self.configuration["pychic_maxLBrownEst"]

        #Seed

        if samples is None:
            samples = 5

        if subset is not None:
            if subset > len(chinput_jiw["baitID"].unique()):
                logger.info("subset > number of baits in data,"+
                            "so used the full dataset.\n")
                subset = None

        if subset == None and samples !=1:
            logger.info("We're using the whole data set to calculate "+
                        "dispersion. There's no reason to sample repeatedly "+
                        "in this case, so overriding brownianNoise.samples to 1.")
            samples = 1

        proxOE = self.readProxOEfile(input_files["poe"])


        #Run this the same as the number of Samples
        dispersion_samples = []
        for i in range(samples):
            dispersion_samples.append(
                self.estimateDispersion(
                            chinput_jiw, proxOE, distFunParams
                            )
            )

        logger.info("Getting consensus dispersion estimate...")

        param_dispersion = np.mean(dispersion_samples)

        chinput_jiw = self.estimateBMean(chinput_jiw, distFunParams)

        return chinput_jiw, param_dispersion

    def getPvals(self, x, dispersion):
        """
        Get the pvalues for each interaction

        Parameters
        ----------
        chinput_jiwb: DataFrame

        dispersion: dict

        Returns
        -------
        chinput_jiwb_pval
        """
        alpha = dispersion

        logger.info("Calculting p-values...")

        ##p-values:
        ##(gives P(X > x-1) = P(X >= x))
        ##Note that the cases Bmean = 0 and Bmean > 0 are considered separately.

        print(x)
        log_p = []

        robjects.numpy2ri.activate()

        delaporte = importr('Delaporte')

        #x["Tmean"].astype(float, inplace=True)

        n_r = robjects.numpy2ri.numpy2ri(np.array(x["N"]))
        tmean_r = robjects.numpy2ri.numpy2ri(np.array(x["Tmean"]))
        bmean_r = robjects.numpy2ri.numpy2ri(np.array(x["Bmean"]))

        robjects.r.assign("N", n_r)
        robjects.r.assign("Tmean", tmean_r)
        robjects.r.assign("Bmean", bmean_r)
        robjects.r.assign("alpha", dispersion)

        robjects.r("""
            Tmean = as.numeric(Tmean)
            """)

        robjects.r("log.p1 <- pdelap(N - 1L, alpha, beta=Bmean/alpha, lambda=Tmean, lower.tail=FALSE, log.p=TRUE)")
        log_p1 = robjects.r("log.p1")
        x["log_p1"] = log_p1


        robjects.r("log.p2 <- ppois(N - 1L, lambda=Tmean, lower.tail=FALSE, log.p=TRUE)")
        log_p2 = robjects.r("log.p2")
        x["log_p2"] = log_p2

        x["log_p"] = np.where(x['Bmean'] < np.finfo(float).eps, x["log_p2"], x["log_p1"])

        x.drop(["log_p1", "log_p2"], axis=1, inplace=True)
        #poisson_sf = np.vectorize(poisson.sf)

        # Large N approximation
        # ---------------------

        ##In rare cases where pdelap returns Infs, estimate the p-value magnitude
        ##using an NB approximation, through method of moments argument
        ##NaNs occur when pdelap() thinks the p-value is negative (since can have 1 - 1 != 0),
        ##thus these are also approximated.

        sel = x[(x["log_p"] == np.nan) | (x["log_p"] == np.inf) | (x["log_p"] == -np.inf)]

        if sel.shape[0] > 0:
            logger.info("Approximating "+ str(len(sel))+ " very small p-values.")


            sel["gamma"] = dispersion * (1+sel["Tmean"]/sel["Bmean"])**2

            sel["gamma"].replace([np.inf, -np.inf], 1e10, inplace=True)

            sel["Bmean+Tmean"] = sel["Bmean"]+ sel["Tmean"]

            gamma_r = robjects.numpy2ri.numpy2ri(np.array(sel["gamma"]))
            bmean_tmean = robjects.numpy2ri.numpy2ri(np.array(sel["Bmean+Tmean"]))
            n_r = robjects.numpy2ri.numpy2ri(np.array(sel["N"]))

            robjects.r.assign("gamma", gamma_r)
            robjects.r.assign("BmeanTmean", bmean_tmean)
            robjects.r.assign("N", n_r)

            robjects.r("log_p <- pnbinom(N - 1L, size=gamma, mu=BmeanTmean, lower.tail=FALSE, log.p=TRUE)")

            log_p = robjects.r("log_p")

            sel["log_p"] = log_p

            sel.sort_values(by=["log_p"])
            x.sort_values(by=["log_p"])

            x.update(sel)

        if x["log_p"].isna().any():
            logger.info("Some log-p-values were NA.")

        return x

    def getAvgFragLength(self, x, rmap, excludeMT=True):
        """
        Get the average fragment

        Parameters
        ----------
        x: dataFrame

        Returns
        -------
        avgFragLen: float
        """
        if excludeMT:
            if "MT" in rmap["chr"].astype(str):
                #Not sure about this MT name
                rmap = rmap[rmap["chr"] != "MT"]

            if "chrMT" in rmap["chr"].astype(str):
                rmap = rmap[rmap["chr"] != "chrMT"]

            chrMAX = rmap.groupby(["chr"], as_index=False)["end"].max()

        avgFragLen = chrMAX["end"].sum()/rmap.shape[0]

        return avgFragLen, rmap


    def getNoOfHypotheses(self, x, rmap, baitmap, avgFragLen):
        """
        Parameter
        ---------
        x: DataFrame

        Returns
        -------
        """
        chrMAX = rmap.groupby(["chr"], as_index=False)["end"].max()

        chromo = chrMAX["chr"]

        ##count # hypothesis
        Nhyp = baitmap.shape[0] * ((2*rmap.shape[0]) - baitmap.shape[0] - 1)/2

        return Nhyp

    def expit(self, x):
        """
        inverse logit function

        Parameter
        ---------
        x : float/int

        Returns
        -------
        exp(x)/(1+exp(x))
        """

        return 1/(1+math.exp(-x))


    def getEtaBar(self, x, rmap, baitmap, avgFragLen):
        """
        get parameters for the experiments

        Parameters
        ----------
        x: Dataframe

        Returns
        -------
        """
        ##1. Collect parameters
        alpha = self.configuration["pychic_weightAlpha"]
        beta = self.configuration["pychic_weightBeta"]
        gamma = self.configuration["pychic_weightGamma"]
        delta = self.configuration["pychic_weightDelta"]

        ##2. Get genomic/fragment map information
        chrMAX = rmap.groupby(["chr"], as_index=False)["end"].max()

        if "MT" in baitmap["chr"].astype(str):
            baitamp = baitmap[baitmap["chr"] != "MT"]

        if "chrMT" in baitmap["chr"].astype(str):
            baitmap = baitmap[baitmap["chr"] != "chrMT"]

        Nhyp = self.getNoOfHypotheses(x, rmap, baitmap, avgFragLen)

        ##3. Calculate eta.bar
        ##Loop, summing contributions of eta

        eta_sigma = 0

        expit = np.vectorize(self.expit)

        for c in rmap["chr"].unique():
            #length of chromosome
            d_c = chrMAX[chrMAX["chr"] == c]
            d_c = d_c["end"]

            nBaits = baitmap[baitmap["chr"] == c]
            n_c = nBaits["chr"].value_counts()

            n_c = int("".join([str(i) for i in n_c]))

            for i in range(1, n_c+1):
                d = round(d_c*i/n_c, 1)

                d_near = min(float(d), float(d_c-d))

                d_other = np.arange(avgFragLen, max(avgFragLen, d_near), avgFragLen)

                d_other2 = np.arange(d_near, d_c-d_near, avgFragLen)

                try:
                    eta_sigma = eta_sigma + 2*sum(expit(alpha + beta*np.log(d_other))) + \
                                sum(expit(alpha + beta*np.log(d_other2)))
                except ValueError:
                    logger.info("Empty array")

        eta_bar = eta_sigma/Nhyp

        return eta_bar

    def getWeights(self, dist, alpha, beta, gamma, delta, eta_bar, expit):
        """
        Calculate the weights
        """

        eta = expit(alpha + beta*np.log(dist))

        a = np.log((self.expit(delta) - self.expit(gamma))*eta + self.expit(gamma))

        b = np.log((self.expit(delta) - self.expit(gamma))*eta_bar) + self.expit(gamma)

        #log_w = np.log(self.expit(delta) - self.expit(gamma))*eta + \
        #               self.expit(gamma) - \
        #        np.log((self.expit(delta) - self.expit(gamma))*eta_bar + self.expit(gamma))

        return a - b

    def getScores(self, x, rmap, baitmap):
        """
        Get the scores from the pvalues

        Parameters
        ----------
        x: DataFrame

        Returns
        -------
        chinput_jiwb_scores: DataFrame
        """
        Set = self.configuration
        avgFragLen, rmap = self.getAvgFragLength(x, rmap)

        eta_bar = self.getEtaBar(x, rmap, baitmap, avgFragLen)

        #Gets weights
        logger.info("Calculating p-values weights...")

        alpha = self.configuration["pychic_weightAlpha"]
        beta = self.configuration["pychic_weightBeta"]
        gamma = self.configuration["pychic_weightGamma"]
        delta = self.configuration["pychic_weightDelta"]

        expit = np.vectorize(self.expit)

        getWeights = np.vectorize(self.getWeights)


        x["log_w"] = getWeights(x["distSign"].abs().replace(np.nan, np.inf),
                                alpha,
                                beta,
                                gamma,
                                delta,
                                eta_bar,
                                expit)

        x["log_q"] = x["log_p"] - x["log_w"]

        logger.info("Calculating scores..")

        ##get score (more interpretable than log.q)
        minval = getWeights(0,
                            alpha,
                            beta,
                            gamma,
                            delta,
                            eta_bar,
                            expit)

        x["score"] = -x["log_q"] - minval

        x["score"] = np.where(x["score"] > 0, x["score"], 0)

        return x


    def print_params(self, pychic_outprefix):
        """
        Print to a file all the parameters used in the experiment

        Parameres
        ---------
        pychic_outprefix: str

        Returns
        -------
        Bool
        """

        file = self.configuration["execution"]+"/"+ \
                   pychic_outprefix+"_params.txt"

        with open(file, "w") as file_out:
            for parameter in self.configuration:
                file_out.write("{}\t{}\n".format(parameter, self.configuration[parameter]))

        return True

    def save_rda(self, x):
        """
        Convert pandas dataframe to an R object and save as a .Rda image

        Parameters
        ----------
        x: DataFrame

        Returns
        -------
        bool
        """
        logger.info("Saving pandas DataFrame as a .Rda object")

        if os.path.isdir(self.configuration["pychic_outprefix"]) == False:
            os.mkdir(self.configuration["pychic_outprefix"])

        import rpy2
        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()

        # write pandas dataframe to an .RData file
        #def save_rdata_file(df, filename):
        r_data = pandas2ri.py2ri(x)
        robjects.r.assign("my_df", r_data)
        robjects.r("save(my_df, file='{}')".format(
            self.configuration["pychic_outprefix"]+"/"+
            self.configuration["pychic_outprefix"]+".Rda"))

        return True

    def cut2(self, num_list, m, onlycuts=True):
        """
        This function si the fixed version of cut2 form R.
        Given a list of numbers and a parameter m, it will bin the list
        having at least m numbers in each bin. The funciton wont guarantee
        that it will be exactly m number of elements in the last bin.
        """
        num_list = sorted(num_list)
        unique_list = list(set(num_list))
        unique_counts = {}
        bins = []

        for number in unique_list:
            unique_counts[number] = num_list.count(number)

        counter = 0

        bins.append(num_list[0])
        for number in unique_list:
            if unique_counts[number] >= m:
                if number in bins:
                    continue
                bins.append(number)
                counter = 0
            else:
                counter += unique_counts[number]
                if counter >= m:
                    bins.append(number)
                    counter = 0
                else:
                    continue

        bins[-1] = num_list[-1]

        return bins

    def exportResults(self, x, outprefix, cutoff, export_format,
                      order, rmap, baitmap):
        """
        print the results in the correct format and order

        Parameter
        ---------
        x: DataFrame
        outprefix: str
        cutoff: int
        format: list
        order: str

        Return
        ------
        bool
        """

        logger.info("Rading the rmap file")

        rmap.columns = ["rChr", "rStart", "rEnd", "otherEndID"]

        baitmap.columns = ["baitChr", "baitStart", "baitEnd",
                           "baitID", "promID"]

        logger.info("Preparing the output table...")

        x = x[x["score"] >= cutoff]

        x = x[["baitID", "otherEndID", "N", "score"]]

        x = pd.merge(x, rmap, how="left", on="otherEndID")

        x = pd.merge(x, baitmap, how="left", on="baitID")

        # note that baitmapGeneIDcol has been renamed into "promID" above
        bm2 = baitmap[["baitID", "promID"]]
        bm2.rename({"promID":"promID_y", "baitID": "otherEndID"}, axis="columns", inplace=True)

        out = pd.merge(x, bm2, how="left", on="otherEndID")

        out["promID_y"] = np.where(out["promID_y"].isna(), ".", out["promID_y"])

        out.rename({"promID":"promID_x"}, axis="columns", inplace=True)

        out = out[["baitChr", "baitStart", "baitEnd", "promID_x",
                   "rChr", "rStart", "rEnd", "otherEndID", "score",
                   "N", "promID_y"]]

        out.columns = ["bait_chr", "bait_start", "bait_end", "bait_name",
                       "otherEnd_chr", "otherEnd_start", "otherEnd_end",
                       "otherEnd_ID", "score", "N_reads", "otherEnd_name"]

        out["N_reads"].fillna(0, inplace=True)
        out["score"] = out["score"].round(2)

        if order == "position":
            out.sort_values(["bait_chr",
                             "bait_start",
                             "bait_end",
                             "otherEnd_chr",
                             "otherEnd_start",
                             "otherEnd_end"],
                            inplace=True)

        elif order == "score":
            out.sort_values(["score"], inplace=True, ascending=False)

        #Modify so it would look at the mt form the begining
        out = out[out["bait_chr"].astype(str).str.lower() != "chrmt"]

        out0 = out
        if "seqMonk" in export_format:
            logger.info("Writing out for seqMonk...")

            out["bait_name"] = out["bait_name"].replace(",","|")

            out["newLineOEChr"] = '\n'+out["otherEnd_chr"].astype(str)

            #out["newLineOEChr"] = out["newLineOEChr"].astype(int)

            out = out[["bait_chr", "bait_start", "bait_end", "bait_name",
                       "N_reads", "score", "newLineOEChr", "otherEnd_start",
                       "otherEnd_end", "otherEnd_name", "N_reads", "score"]]

            #print(out)
            import csv
            out.to_csv(self.configuration["pychic_outprefix"]+"_seqmonk.txt",
                       sep="\t",
                       header=False,
                       index=False,
                       doublequote=False,
                       quoting=csv.QUOTE_NONE,
                       quotechar="", escapechar="\t"
                       )

        if "interBed" in export_format:
            logger.info("writting out interBed...")
            out = out0[["bait_chr", "bait_start", "bait_end", "bait_name",
                        "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                        "N_reads", "score"]]

            out.to_csv(self.configuration["pychic_outprefix"]+".ibed",
                       sep="\t",
                       index=False
                       )

        if "washU_text" in export_format or "washU_track" in export_format:
            logger.info("Preprocessing for WashU outputs...")

            out = out0[["bait_chr", "bait_start", "bait_end",
                        "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                        "score"]]

            if str(out.iloc[0:1, 0][:3]) != "chr":
                out["bait_chr"] = "chr"+out["bait_chr"].astype(str)
            if str(out.iloc[0:1, 3][:3]) != "chr":
                out["otherEnd_chr"] = "chr"+out["otherEnd_chr"].astype(str)

            # Bait to bait interactions can be asymmetric in terms of score.
            # Here, we find asymmetric interactions and delete the minimum score

            x.sort_values("score").drop_duplicates(subset=["baitID","otherEndID"], keep="last")

        if "washU_text" in export_format:
            logger.info("Writing out text file for WashU browser upload...")
            if os.path.isdir(self.configuration["pychic_outprefix"]) == False:
                os.mkdir(self.configuration["pychic_outprefix"])

            res = pd.DataFrame()
            res["1bait"] = out["bait_chr"].astype(str)+","+out["bait_start"].astype(str)+","+out["bait_end"].astype(str)
            res["2bait"] = out["otherEnd_chr"].astype(str)+","+out["otherEnd_start"].astype(str)+","+out["otherEnd_end"].astype(str)
            res["score"] = out["score"]

            name = self.configuration["pychic_outprefix"]+"/"+ \
                        self.configuration["pychic_outprefix"]+"_washU_text.txt"

            res.to_csv(name, sep="\t", header=False, index=False)
        """
        if "washU_track" in export_format:
            logger.info("Writing out track for WashU browser...")

            ## this format requires a duplicate of each row, with bait/otherEnd reversed
            out["i"] = np.ones(len(out))
            appendOut =  out[["otherEnd_chr", "otherEnd_start", "otherEnd_end", "bait_chr", "bait_start", "bait_end", "otherEnd_name", "score", "i"]]
            appendOut.columns = out.columns

            out = out.append(appendOut)
            out.sort_values(["bait_chr", "bait_start"], inplace=True)

            print(out)
        """

    def plotBaits(self, x, dispersion):
        """
        This function generates a pdf with 10 random baits
        """
        baits = x["baitID"].unique()
        sns.set_style("ticks")

        if len(baits) < 6:
            baits = baits
            row = int(len(baits)/2)
            col = 2
        else:
            baits = sample(list(baits), 6)
            row = 3
            col = 2

        x = x[x["baitID"].isin(baits)]

        baits = list(baits)

        fig, axs = plt.subplots(row, col, figsize=(15,15))

        index = 0
        for i in range(col):
            for j in range(row):

                x_temp = x[x["baitID"] == baits[index]]
                x_temp = x_temp[x_temp["distSign"].abs() < 800000]
                x_temp.sort_values(by="distSign", inplace=True)

                axs[j, i].set_ylim(0, max(x_temp["N"]+1))
                axs[j, i].set_xlim(-800000, 800000)

                x_red = x_temp[x_temp["score"] >= 5]
                x_blue = x_temp[(x_temp["score"] >= 3) & (x_temp["score"] < 5)]
                x_black = x_temp[x_temp["score"] < 3]

                axs[j, i].scatter(x_red["distSign"], x_red["N"], s=15, c="red", label= "Significant interactions")
                axs[j, i].scatter(x_blue["distSign"], x_blue["N"], s=15, c="blue", label = "Sub-threshold interactions")
                axs[j, i].scatter(x_black["distSign"], x_black["N"], s=15,  facecolors='none', edgecolors='black', label="Non significant interactions")

                axs[j, i].plot([0, max(x_temp["N"])], c="black")
                axs[j, i].set_xlabel("Genomic distance from viewpoint")
                axs[j, i].set_ylabel("Number of interactions")

                significant = x_temp["Bmean"].astype(float)+1.96*np.sqrt(
                    x_temp["Bmean"].astype(float)+x_temp["Bmean"].astype(float)**2/dispersion)

                axs[j, i].fill_between(x_temp["distSign"],0,
                                       significant,
                                       facecolor='grey', alpha=0.5, label="Background model")

                axs[j, i].legend(loc='upper left')

                index += 1

        plt.show()
        #plt.savefig(self.configuration["pychic_outprefix"]+"_examples.pdf",
        #            quality=95, dpi="figure", facecolor='w', edgecolor='w',
        #            papertype=None,
        #            box_inches='tight', pad_inches="tight")


    def run(self, input_files, metadata, output_files):
        """
        Function that runs and pass the parameters to PyCHiC

        Parameters
        ----------
        input_files : dict
            pychic_rmap
            pychic_baitmap
            pychic_poe
            pychic_npbp
            pychic_npb
            pychic_settings-file
            pychic_chinput

        metadata : dict
            pychic_cutoff
            pychic_export_format
            pychic_export_order
            pychic_examples_prox_dist
            pychic_examples_full_range

        output_files : dict
            pychic_output_dir: path to the output folder
            pychic_outprefix: output prefix

        Returns
        -------
        output_files : dict
        output_metadata : dict
        """

        if os.path.isdir(self.configuration["pychic_output_dir"]) is False:
            os.makedirs(self.configuration["pychic_output_dir"])

        for test_file in input_files:
            if self.exist_file(input_files[test_file], test_file) is False:
                return False

        self.checks(
            self.configuration["pychic_export_format"],
            self.configuration["pychic_order"],
            input_files["RMAP"],
            input_files["BAITMAP"],
            input_files["nbpb"],
            input_files["npb"],
            input_files["poe"],
            input_files["pychic_settings_file"]
        )


        if len(input_files["chinput"]) == 1:
            chinput_filtered = self.readSample(input_files["chinput"],
                                               self.configuration["pychic_bam"],
                                               input_files["RMAP"],
                                               input_files["BAITMAP"])

            #chinput_filtered.to_csv("chinput_filteres.csv", sep="\t", index=False)

            #npb = self.prepare_design(input_files["npb"], "baitID")

        else:
            chinputs_filtered = {}
            for i in range(len(input_files["chinput"])):
                new_chinput = self.readSample(input_files["chinput"][i],
                                              self.configuration["pychic_bam"],
                                              input_files["RMAP"],
                                              input_files["BAITMAP"])

                chinputs_filtered[str(i)] = new_chinput

            chinput_filtered, npb = self.merge_chinputs(chinputs_filtered, input_files["npb"])

        logger.info("\nRunning normaliseBaits")

        chinput_j = self.normaliseBaits(chinput_filtered, \
                                       input_files["npb"])


        chinput_ji = self.normaliseOtherEnds(chinput_j,
                                             input_files["nbpb"]
                                            )


        logger.info("\n Running estimateTechicalNoise")

        chinput_jiw = self.estimateTechnicalNoise(chinput_ji,
                                                  input_files["RMAP"],
                                                  input_files["BAITMAP"])


        distFunParams = self.estimateDistFun(chinput_jiw)

        chinput_jiwb, dispersion = self.estimateBrownianComponent(chinput_jiw, distFunParams)

        chinput_jiwb_pval = self.getPvals(chinput_jiwb, dispersion)

        rmap_df = pd.read_csv(input_files["RMAP"],
                              sep="\t",
                              names=["chr", "start", "end", "ID"])

        baitmap_df = pd.read_csv(input_files["BAITMAP"],
                                 header=None,
                                 sep="\t")

        if len(baitmap_df.columns) == 4:
            baitmap_df.columns = ["chr", "start", "end", "ID"]
        else:
            baitmap_df.columns = ["chr", "start", "end", "ID", "feature"]


        chinput_jiwb_scores = self.getScores(chinput_jiwb_pval,
                                             rmap_df,
                                             baitmap_df
                                            )

        self.print_params(self.configuration["pychic_outprefix"])

        if self.configuration["pychic_Rda"] == "True":
            self.save_rda(chinput_jiwb_scores)

        #chinput_jiwb_scores = pd.read_csv("chinput_jiwb_scores.txt", sep="\t")

        #dispersion = 2.5563913

        #self.plotBaits(chinput_jiwb_scores, dispersion)

        self.exportResults(chinput_jiwb_scores,
                           self.configuration["pychic_outprefix"],
                           self.configuration["pychic_cutoff"],
                           self.configuration["pychic_export_format"],
                           self.configuration["pychic_order"],
                           rmap_df,
                           baitmap_df
                          )

if __name__ == "__main__":

    path = "../../tests/data/test_run_chicago/data_chicago/"

    input_files = {
         "RMAP" : path +"h19_chr20and21.rmap",
         "BAITMAP" : path +"h19_chr20and21.baitmap",
         "nbpb" : path +"h19_chr20and21.nbpb",
         "npb" : path +"h19_chr20and21.npb",
         "poe" : path +"h19_chr20and21.poe",
         "pychic_settings_file" : path +"sGM12878.settingsFile",
         "chinput" :  [
                      path + "GM_rep1.chinput",
                      #path + "GM_rep2.chinput",
                      #path + "GM_rep3.chinput"
                      ]
    }

    metadata = {}

    output_files = {

    }

    configuration = {
        "pychic_cutoff" : 5,
        "pychic_Rda" : "False",
        "execution" : ".",
        "pychic_outprefix" : "out_test",
        "pychic_export_format" : ["washU_text"],
        "pychic_order" : "score",
        "pychic_output_dir" : path +"output_pyCHiC",
        "pychic_rmapfile" : "NA",
        "pychic_baitmapfile" :"NA",
        "pychic_nperbinfile" : "NA",
        "pychic_nbaitsperbinfile" : "NA",
        "pychic_proxOEfile" : "NA",
        "pychic_Ncol" : "N",
        "pychic_baitmapFragIDcol" : 4,
        "pychic_baitmapGeneIDcol" : 5,
        "pychic_maxLBrownEst" : 1500000.0,
        "pychic_minFragLen" : 150, # minimun OE fragment lenght in bps
        "pychic_maxFragLen" : 40000, # maximun OE fragment lenght in bps
        "pychic_minNPerBait" : 250,
        "pychic_binsize" : 20000,
        "pychic_removeAdjacent" : True,
        "pychic_adjBait2bait" : True,
        "pychic_tlb_filterTopPercent" : 0.01,
        "pychic_tlb_minProxOEPerBin" : 150,
        "pychic_tlb_minProxB2BPerBin" : 15,
        "pychic_techNoise_minBaitsPerBin" : 150,
        "pychic_brownianNoise_samples" : 5,
        "pychic_brownianNoise_subset" : 500,
        "pychic_brownianNoise_seed" : "NA",
        "pychic_baitIDcol" : "baitID",
        "pychic_otherEndIDcol" : "otherEndID",
        "pychic_otherEndLencol" : "otherEndLen",
        "pychic_distcol" : "distSign",
        "pychic_weightAlpha" : 34.1157346557331,
        "pychic_weightBeta" : -2.58688050486759,
        "pychic_weightGamma" : -17.1347845819659,
        "pychic_weightDelta" : -7.07609245521541,
        "pychic_bam" : "sample_549_lane7_NoIndex_L007_R1_2.hicup.captured"
    }

    pyCHiC_obj = pyCHiC(configuration)
    pyCHiC_obj.run(input_files, metadata, output_files)


"""
Write to output the parameters
"""