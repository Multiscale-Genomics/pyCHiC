
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

#from pytadbit.parsers.hic_parser import load_hic_data_from_reads
#from pytadbit.parsers.map_parser import parse_map
#import pickle
from pytadbit.utils.file_handling import mkdir, which
from collections                  import OrderedDict
#from pytadbit.mapping.filter      import MASKED
from distutils.version            import LooseVersion
from subprocess                   import Popen, PIPE
import sys

MASKED = {1 : {'name': 'self-circle'       , 'reads': 0},
          2 : {'name': 'dangling-end'      , 'reads': 0},
          3 : {'name': 'error'             , 'reads': 0},
          4 : {'name': 'extra dangling-end', 'reads': 0},
          5 : {'name': 'too close from RES', 'reads': 0},
          6 : {'name': 'too short'         , 'reads': 0},
          7 : {'name': 'too large'         , 'reads': 0},
          8 : {'name': 'over-represented'  , 'reads': 0},
          9 : {'name': 'duplicated'        , 'reads': 0},
          10: {'name': 'random breaks'     , 'reads': 0},
          11: {'name': 'trans-chromosomic' , 'reads': 0}}



def _map2sam_chicago(line, flag=0):
    """
    translate map + flag into hic-sam (two lines per contact)
    only loses RE sites, that can be added back later
    63% of the size using RE sites, and 72% of the generation time
    """
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

def bed2D_to_BAMhic(infile, valid, ncpus, outbam, frmt ='chicago', masked=None, samtools='samtools'):
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
       - pseudo CIGAR with sequence length and info about current copy (P: first copy, S: second copy)
       - chromosome ID of the second pair of the contact
       - genomic position of the second pair of the contact
       - mapped length of the second pair of the contact
       - sequence is missing (*)
       - quality is missing (*)
       - TC tag indicating single (1) or multi contact (3 6 ... number being the number of times a given sequenced fragment is involved in a pairwise contact)
       - S1 and S2 tags are the strand orientation of the left and right read-end
    Each pair of contacts produces two lines in the output BAM
    """
    samtools = which(samtools)
    if not samtools:
        raise Exception('ERROR: samtools is needed to save a compressed '
                        'version of the results. Check '
                        'http://samtools.sourceforge.net/ \n')

    # define filter codes
    filter_keys = OrderedDict()
    for k in MASKED:
        filter_keys[MASKED[k]['name'].replace(' ', '-')] = 2 ** (k - 1)

    output = ''

    # write header
    output += ("\t".join(("@HD" ,"VN:1.0", "SO:coordinate")) + '\n')
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
    output += ("\t".join(("@CO" ,"XA:i", "Number of time a sequenced fragment is involved in a pairwise contact\n")))
    output += ("\t".join(("@CO" ,"NM:i", " Edit distance to the reference, including ambiguous bases but ",
                          "excluding clipping\n")))
    output += ("\t".join(("@CO" ,"MD:Z", " String for mismatching positions\n")))
    output += ("\t".join(("@CO" ,("Each read is duplicated: once starting with the "
                                  "left read-end, once with the right read-end\n"))))
    output += ("\t".join(("@CO" , (" the order of RE sites and strands changes consequently "
                                   "depending on which read-end comes first ("
                                   "when right end is first: E3 E4 E1 E2)\n"))))
    output += ("\t".join(("@CO" ,(" CIGAR code contains the length of the "
                                  "1st read-end mapped and 'P' or 'S' "
                                  "if the copy is the first or the second\n"))))
    output += ("\t".join(("@CO" ,"E1:i", "Position of the left RE site of 1st read-end\n")))
    output += ("\t".join(("@CO" ,"E2:i", "Position of the right RE site of 1st read-end\n")))
    output += ("\t".join(("@CO" ,"E3:i", "Position of the left RE site of 2nd read-end\n")))
    output += ("\t".join(("@CO" ,"E4:i", "Position of the right RE site of 2nd read-end\n")))
    output += ("\t".join(("@CO" ,"S1:i", "Strand of the 1st read-end (1: positive, 0: negative)\n")))
    output += ("\t".join(("@CO" ,"S2:i", "Strand of the 2nd read-end  (1: positive, 0: negative)\n")))
    print(1)
    # open and init filter files
    if not valid:
        filter_line, filter_handler = get_filters(infile, masked)
    fhandler.seek(pos_fh)
    # check samtools version number and modify command line
    version = LooseVersion([l.split()[1]
                            for l in Popen(samtools, stderr=PIPE).communicate()[1].split('\n')
                            if 'Version' in l][0])
    version = "1.9"
    pre = '-o' if version >= LooseVersion('1.3') else ''
    print(2)
    proc = Popen(samtools + ' view -Shb -@ %d - | samtools sort -@ %d - %s %s' % (
        ncpus, ncpus, pre,
        outbam + '.bam' if  version >= LooseVersion('1.3') else ''),  # in new version '.bam' is no longer added
                 shell=True, stdin=PIPE)
    proc.stdin.write(output)
    print(3)
    if frmt == 'chicago':
        map2sam = _map2sam_chicago
    elif frmt == 'mid':
        map2sam = _map2sam_mid
    elif frmt == 'long':
        map2sam = _map2sam_long
    else:
        map2sam = _map2sam_short

    if valid:
        for line in fhandler:
            flag = 0
            # get output in sam format
            proc.stdin.write(map2sam(line, flag))
        print(4)
    else:
        for line in fhandler:
            flag = 0
            # check if read matches any filter
            rid = line.split("\t")[0]
            for i in filter_line:
                if filter_line[i] == rid:
                    flag += filter_keys[i]
                    try:
                        filter_line[i] = filter_handler[i].next().strip()
                    except StopIteration:
                        pass
            # get output in sam format
            proc.stdin.write(map2sam(line, flag))
    print(5)
    proc.stdin.close()
    proc.wait()

    # Index BAM
    #_ = Popen(samtools + ' index %s.bam' % (outbam), shell=True).communicate()
    print(6)
    # close file handlers
    fhandler.close()
    if not valid:
        for i in filter_handler:
            filter_handler[i].close()
    print(7)

if "__main__" == __name__:
    print(0)
    print(sys.argv)
    bed2D_to_BAMhic(sys.argv[1],"valid",int(sys.argv[2]),sys.argv[3])

    #bed2D_to_BAMhic("/home/pablo/MuG/data/trim_galore_output/tadbit_truncated/03_filtered_reads/valid_r1-r2_intersection_b51cdf1282.tsv",
                    #     "valid", 2, "/home/pablo/MuG/data/trim_galore_output/tadbit_truncated/")


#REMEMBER TO CALL SAMTOOLS TO SORT THE OUTPUT
