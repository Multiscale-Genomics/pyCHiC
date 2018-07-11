
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
import subprocess
from utils import logger
from shutil import copy
import shlex

from pytadbit.parsers.hic_parser import load_hic_data_from_reads
from pytadbit.parsers.map_parser import parse_map
import pickle
from pytadbit.utils.file_handling import mkdir, which
from collections                  import OrderedDict
from subprocess                   import Popen, PIPE


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

class bed2bam(Tool):
    """
    This class contain functions to convert a bed file bam_out from
    fatq2bed.py to bam file compatible with CHiCAGO
    """
    def __init__(self, configuration=None):
        """
        Initiate the tool

        Parameters
        ----------
        configuration: dict
         contain info to run the functions of the class
        """

        logger.info("Initiating bed2chicago")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)


    def _map2sam_chicago(self, line, flag=0):
        """
        translate map + flag into hic-sam (two lines per contact)
        only loses RE sites, that can be added back later
        63% of the size using RE sites, and 72% of the generation time
        """
        logger.info("_map2sam_chicago from bed2bam function running")

        (qname,
         rname, pos, s1, l1, _, _,
         rnext, pnext, s2, l2, _) = line.strip().split('\t', 11)

        # multicontact?
        try:
            tc = qname.split('#')[1].split('/')[1]
        except IndexError:
            tc = '1'

        # Make flag
        # First add read paired and read mapped proper in pair
        flag = 1 + 2
        # Add info about strand
        s1 = (int(s1) - 1)
        s2 = (int(s2) - 1)
        flag1 = flag + (-16 * s1) + (-32 * s2)
        flag2 = flag + (-16 * s2) + (-32 * s1)
        # Then we add pair info
        flag1 += 64
        flag2 += 128


        r1r2 = ('{0}\t{1}\t{2}\t{3}\t255\t{4}M\t{6}\t{7}\t0\t*\t*\t'
                'XA:i:{8}\tMD:Z:{4}\tNM:i:0\n'
                '{0}\t{9}\t{6}\t{7}\t255\t{5}M\t{2}\t{3}\t0\t*\t*\t'
                'XA:i:{8}\tMD:Z:{5}\tNM:i:0\n').format(
                    qname,               # 0
                    flag1,               # 1
                    rname,               # 2
                    pos,                 # 3
                    l1,                  # 4
                    l2,                  # 5
                    rnext,               # 6
                    pnext,               # 7
                    tc,                  # 8
                    flag2)               # 9
        return r1r2

    @task(returns=bool, infile=FILE_IN, valid=IN, ncpus=IN,
          outbam=FILE_OUT)
    def bed2D_to_BAMhic(self, infile, valid, ncpus, outbam):

        """
        function adapted from Enrique Vidal <enrique.vidal@crg.eu> scipt to convert
        2D beds into compressed BAM format.
        Gets the *_both_filled_map.tsv contacts from TADbit (and the corresponding
        filter files) and outputs a modified indexed BAM with the following fields:
           - read ID
           - filtering flag (see codes in header)
           - chromosome ID of the first pair of the contact
           - genomic position of the first pair of the contact
           - MAPQ set to 0
           - pseudo CIGAR with sequence length and info about current
                copy (P: first copy, S: second copy)
           - chromosome ID of the second pair of the contact
           - genomic position of the second pair of the contact
           - mapped length of the second pair of the contact
           - sequence is missing (*)
           - quality is missing (*)
           - TC tag indicating single (1) or multi contact
                (3 6 ... number being the number of times a given sequenced
                fragment is involved in a pairwise contact)
           - S1 and S2 tags are the strand orientation of the left and right read-end
        Each pair of contacts produces two lines in the output BAM
        """
        frmt = 'chicago'
        samtools = 'samtools'

        try:
            logger.info("bed2D_to_BAMhic from bed2bam function running")

            masked = {1 : {'name': 'self-circle', 'reads': 0},
                      2 : {'name': 'dangling-end', 'reads': 0},
                      3 : {'name': 'error', 'reads': 0},
                      4 : {'name': 'extra dangling-end', 'reads': 0},
                      5 : {'name': 'too close from RES', 'reads': 0},
                      6 : {'name': 'too short', 'reads': 0},
                      7 : {'name': 'too large', 'reads': 0},
                      8 : {'name': 'over-represented', 'reads': 0},
                      9 : {'name': 'duplicated', 'reads': 0},
                      10: {'name': 'random breaks', 'reads': 0},
                      11: {'name': 'trans-chromosomic', 'reads': 0}}

            samtools = which(samtools)
            if not samtools:
                raise Exception('ERROR: samtools is needed to save a compressed '
                                'version of the results. Check '
                                'http://samtools.sourceforge.net/ \n')

            # define filter codes
            filter_keys = OrderedDict()
            for k in masked:
                filter_keys[masked[k]['name'].replace(' ', '-')] = 2 ** (k - 1)

            output = ''

            # write header
            output += ("\t".join(("@HD", "VN:1.0", "SO:coordinate")) + '\n')
            fhandler = open(infile)
            line = fhandler.next()
            # chromosome lengths
            pos_fh = 0

            # In samtools this is sorted alphabetically
            header = {}
            while line.startswith('#'):
                (_, _, cr, ln) = line.replace("\t", " ").strip().split(" ")
                header[cr] = ("\t".join(("@SQ", "SN:" + cr, "LN:" + ln)) + '\n')
                #output += ("\t".join(("@SQ", "SN:" + cr, "LN:" + ln)) + '\n')
                pos_fh += len(line)
                line = fhandler.next()
            hdrOrdr = sorted(header.keys())
            for key in hdrOrdr:
                output += header[key]

            # filter codes
            for i in filter_keys:
                output += ("\t".join(("@CO", "filter:" + i, "flag:" + str(filter_keys[i]))) + '\n')

            # tags
            output += ("\t".join(("@CO", "XA:i",
                                  "Number of time a sequenced fragment"
                                  " is involved in a pairwise contact\n")))
            output += ("\t".join(("@CO", "NM:i",
                                  " Edit distance to the reference, including ambiguous bases but ",
                                  "excluding clipping\n")))
            output += ("\t".join(("@CO", "MD:Z",
                                  " String for mismatching positions\n")))
            output += ("\t".join(("@CO", ("Each read is duplicated: once starting with the "
                                          "left read-end, once with the right read-end\n"))))
            output += ("\t".join(("@CO", (" the order of RE sites and strands changes consequently "
                                          "depending on which read-end comes first ("
                                          "when right end is first: E3 E4 E1 E2)\n"))))
            output += ("\t".join(("@CO", (" CIGAR code contains the length of the "
                                          "1st read-end mapped and 'P' or 'S' "
                                          "if the copy is the first or the second\n"))))
            output += ("\t".join(("@CO", "E1:i",
                                  "Position of the left RE site of 1st read-end\n")))
            output += ("\t".join(("@CO", "E2:i",
                                  "Position of the right RE site of 1st read-end\n")))
            output += ("\t".join(("@CO", "E3:i",
                                  "Position of the left RE site of 2nd read-end\n")))
            output += ("\t".join(("@CO", "E4:i",
                                  "Position of the right RE site of 2nd read-end\n")))
            output += ("\t".join(("@CO", "S1:i",
                                  "Strand of the 1st read-end (1: positive, 0: negative)\n")))
            output += ("\t".join(("@CO", "S2:i",
                                  "Strand of the 2nd read-end  (1: positive, 0: negative)\n")))

            fhandler.seek(pos_fh)

            try:
                cmd_1 = samtools + " view -Shb -@ {} -".format(ncpus)
                cmd_2 = samtools + " sort -@ {} - {} {}".format(ncpus, "-o", outbam+".tmp")

                p1 = Popen(shlex.split(cmd_1), stdin=PIPE, stdout=PIPE)
                p2 = Popen(shlex.split(cmd_2), stdin=p1.stdout)
                p1.stdin.write(output)
            except IOError:
                logger.fatal("Popen function does not work for samtools =( ")
                return False

            """
            proc = Popen(samtools + ' view -Shb -@ %d - | samtools sort -@ %d - %s %s' % (
                         ncpus, ncpus, "-o", outbam+".tmp"),  # in new version '.bam' is no longer added
                         shell=True, stdin=PIPE)
            proc.stdin.write(output)
            """
            if frmt == 'chicago':
                map2sam = self._map2sam_chicago

            if valid:
                for line in fhandler:
                    flag = 0
                    # get output in sam format
                    p1.stdin.write(map2sam(line, flag))

            p1.stdin.close()
            p1.wait()

            output = p2.communicate()[0]
            # close file handlers
            fhandler.close()

            try:
                with open(outbam+".tmp", "r") as f_in:
                    with open(outbam, "w") as f_out:
                        f_out.write(f_in.read())
                logger.info("tmp bam file converted to bam_out")
                os.remove(outbam+".tmp")
                return True

            except IOError:
                logger.fatal("temporary file not converted to output bam file")
                return False

        except IOError:
            return False

    @task(returns=bool, bam_out=FILE_IN,
          bam_out_sorted=FILE_OUT)
    def sort_bam_out(self, bam_out, bam_out_sorted):
        """
        This function sort the bam_out using samtools

        Parameters
        ----------
        bam_out: str
            path to bam_out directory and file
        """
        args = ["samtools", "sort",
                "-n", bam_out]

        logger.info("samtools args:"+ " ".join(args))

        try:
            with open(bam_out_sorted, "w") as f_out:
                process = subprocess.Popen(
                    ' '.join(args),
                    shell=True,
                    stdout=f_out, stderr=f_out
                    )
                process.wait()
            return True

        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}): {1}\n{2}".format(
                msg.errno, msg.strerror, args))
            return False

    def run(self, input_files, metadata, output_files):
        """
        This function runs the wrapper_bed2chicago function
        and produce the bam_out

        Parameters
        ----------
        input_files: dict
            deb
            ncpus
        metadata: dict
        output_files: dict
            bam_out

        Returns
        -------
        results: bool
        output_metadata:dict
        """
        output_dir = os.path.split(output_files["bam_out"])[0]

        if os.path.isdir(output_dir) is False:
            logger.info("creating output directory")
            os.mkdir(output_dir)

        ncpus = self.configuration["ncpus"]

        #results = self.wrapper_bed2bam(
        #    input_files["bed"],
        #    output_files["bam_out"])
        results = self.bed2D_to_BAMhic(
            input_files["bed"],
            "valid",
            2,
            output_files["bam_out"])

        results = compss_wait_on(results)

        if results is True:
            sorted_results = self.sort_bam_out(
                output_files["bam_out"], output_files["bam_out_sorted"])

            sorted_results = compss_wait_on(sorted_results)

        output_metadata = {
            "bam_out": Metadata(
                data_type="TXT",
                file_type="bam",
                file_path=output_files["bam_out"],
                sources=metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                ),
            "bam_out_sorted": Metadata(
                data_type="TXT",
                file_type="bam",
                file_path=output_files["bam_out_sorted"],
                sources=metadata["bed"].file_path,
                taxon_id=9606,
                meta_data=""
                ),
        }

        return output_files, output_metadata
