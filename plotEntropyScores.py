#!/lab/64/bin/python
# plotEntropyScores.py 
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Will plot a histogram of the Shannon Entropy scores for all junctions (to do).

It will also return a file containing the entropy score for all junctions and
the offset positions for every read along with the entropy.
"""

import sys
import optparse 
import math

import pysam
from preProcess_getASEventReadCounts import JcnInfo, convert2SAMLine, getForcedJunctions

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy
#############
# CONSTANTS #
#############

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("-s",
                          dest="sam_file",
                          type="string",
                          help="BAM/SAM file which contains junction reads.",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_name",
                          type="string",
                          help="Prefix name for output files.",
                          default=None)
    opt_parser.add_option("--known_junctions",
                          dest="known_junctions",
                          type="string",
                          help="""File containing intron coordinates of known/annotated
                                  introns. If this is given, the entropy scores
                                  will be flagged as (K)nown or (N)ovel
                                  junctions.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-s")
    opt_parser.check_required("-o")

    isBam = False
    if options.sam_file.endswith(".bam"):
        sam_file = pysam.Samfile(options.sam_file, "rb")
        isBam = True
    else:
        sam_file = open(options.sam_file)
    output_name = options.output_name

    known_junctions = None
    if options.known_junctions:
        known_junctions = getForcedJunctions(options.known_junctions)
    

    jcn2JcnInfoDict, jcn2type = parseSAMFile(sam_file, known_junctions, isBam)

#    all_entropies = []

    output_file = open(output_name + "_entropy_scores.txt", "w")
    offset_file = open(output_name + "_entropy_offset.txt", "w")

    entropy_scores = []    
    if known_junctions:
        novel_entropy_scores = []

    for jcn in jcn2JcnInfoDict:
        pos2count = getPos2Count(jcn2JcnInfoDict[jcn].block_list)
        totalCount = len(jcn2JcnInfoDict[jcn].block_list)

        entropy = getShannonIndex(pos2count, totalCount)

#        all_entropies.append(entropy)
        if known_junctions:
            output_file.write("%s\t%.3f\t%s\n" % (jcn, entropy, jcn2type[jcn]))
            if jcn2type[jcn] == "N":
                novel_entropy_scores.append(jcn2type[jcn])
            else:
                entropy_scores.append(jcn2type[jcn])
        else:
            output_file.write("%s\t%.3f\n" % (jcn, entropy))

        for upstr_overhang in jcn2JcnInfoDict[jcn].block_list:
            if known_junctions:
                offset_file.write("%.3f\t%d\t%s\n" % (entropy, upstr_overhang, jcn2type[jcn]))
            else:
                offset_file.write("%.3f\t%d\n" % (entropy, upstr_overhang))

    output_file.close()
    offset_file.close()

    # Create entropy score distribution
    fig = plt.figure()
    
    plt.hist(numpy.array(entropy_scores), 
             bins=20,
             range=[0,10],
             histtype='stepfilled',
             normed=True,
             color = 'b',
             alpha = 0.25,
             label="junction entropy")

    if known_junctions:
        plt.hist(numpy.array(novel_entropy_scores), 
                 bins=20,
                 range=[0,10],
                 histtype='stepfilled',
                 normed=True,
                 color = 'r',
                 alpha = 0.25,
                 label="junction entropy")

    plt.xlabel("Entropy Score")
    plt.ylabel("Probability")
    plt.legend()

    fig.savefig("%s_junction_entropy_distribution.png" % output_name)
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getPos2Count(blocklist):
    pos2count = {}

    for block in blocklist:
        if block in pos2count:
            pos2count[block] += 1
        else:
            pos2count[block] = 1

    return pos2count

def getShannonIndex(pos2countDict, totalCount):
    """
    Calculates a Shannon-Wiener Index for the positions and counts.  It treats
    each position like a unique species.
    """

    if totalCount == 0:
        return 0

    summation = 0

    for pos in pos2countDict:
        p = float(pos2countDict[pos]) / totalCount

        if p == 0.0:
            continue

        summation += p * math.log(p, 2)

    if summation == 0:
        return 0

    return -summation

def getType(jcn_str, known_junctions):
    if jcn_str in known_junctions:
        return "K"

    return "N"
    

def parseSAMFile(sam_file, known_junctions, isBam):
    jcn2JcnInfo = {}
    jcn2type = {}

    insertionFlag = False
    deletionFlag = False
    softclipFlag = False
    hardclipFlag = False

    for line in sam_file:
        if isBam:
            line = convert2SAMLine(sam_file, line)

        line = formatLine(line)

        # Ignore headers
        if line.startswith("@"):
            continue

        sam_elems = line.split("\t")

        q_name = sam_elems[0]

        flag = int(sam_elems[1])
        chr = sam_elems[2]

        if not chr.startswith("chr"):
            chr = "chr" + chr

        chr_start = int(sam_elems[3])

        cigar = sam_elems[5]
        tags = sam_elems[11:]

        m_count = cigar.count("M")

        # Check if it is a genome read
        if m_count > 1: # A JUNCTION READ
            n_count = cigar.count("N")

            i_count = cigar.count("I")
            if i_count > 0:
                if not insertionFlag:
                    print "Not supporting insertions, yet e.g., %s" % cigar
                insertionFlag = True
                continue

            d_count = cigar.count("D")
            if d_count > 0:
                if not deletionFlag:
                    print "Not supporting deletions, yet e.g., %s" % cigar
                deletionFlag = True
                continue

            s_count = cigar.count("S")
            if s_count > 0:
                if not softclipFlag:
                    print "Not supporting softclipping, yet e.g., %s" % cigar
                softclipFlag = True
                continue

            h_count = cigar.count("H")
            if h_count > 0:
                if not hardclipFlag:
                    print "Not supporting hardclipping, yet e.g., %s" % cigar
                hardclipFlag = True
                continue


            if n_count == 0:
                print "Expecting a junction read: %s" % cigar
                continue

            n_split = cigar.split("N")

            # Get the downstrm length
            downstr_len = int(n_split.pop().rstrip("M"))

            # A list to hold the information about each intron. Used in cases
            # where a read aligns to multiple junctions.
            introns_info = []

            first_chr_start = chr_start

            # Get first intron information which also will be used for the
            # upstream length
            upstr_len, intron_len = map(int, n_split.pop(0).split("M"))
            introns_info.append((chr_start, intron_len, None))
            # Updating the chr_start
            chr_start = chr_start + upstr_len + intron_len

            # Get remaining information from additional introns if there are
            # any
            # 3rd element is used in calculating entropy
            for remaining_intron in n_split:
                exon_len, intron_len = map(int, remaining_intron.split("M"))
                introns_info.append((chr_start + exon_len - upstr_len,
                                     intron_len,
                                     chr_start - first_chr_start))
                chr_start = chr_start + exon_len + intron_len


            jcn_tag = None
            jcn_strand = None
            type = None
            for tag in tags:
                if tag.startswith("Y0"):
                    jcn_tag = tag[5:]

                if tag.startswith("XS"):
                    tag, almost_strand = tag.split("A")
                    jcn_strand = almost_strand.lstrip(":")

                    if jcn_strand != "+" and jcn_strand != "-" and jcn_strand != ".":
                        print "Error in strand information for: %s" % line
                        sys.exit(1)

            # Now insert all introns into jcn dictionary.
            for intron_info in introns_info:
                chr_start = intron_info[0]
                intron_len = intron_info[1]

                total_len = upstr_len + intron_len + downstr_len

                jcn_str = None
                if not jcn_tag:
                    # Create a junction string based on the 1-based junction
                    # coordinate
                    jcn_str = "%s_%d_%d" % (chr,
                                            chr_start + upstr_len,
                                            chr_start + upstr_len + intron_len - 1)
                else:
                    # Need to make multiple jcn_str for each intron
                    jcn_str = "%s_%s_%d_%d" % (jcn_tag,
                                               chr,
                                              chr_start + upstr_len,
                                              chr_start + upstr_len + intron_len - 1)

                if known_junctions:
                    type = getType(jcn_str, known_junctions)

                if not jcn_strand:
                    jcn_strand = "."

                # Get BED format information 
                chromStart = chr_start - 1
                chromEnd = chromStart + total_len

                # Now add junction to dictionary
                if jcn_str in jcn2JcnInfo:
                    jcn2JcnInfo[jcn_str].updateJcnInfo(jcn_str,
                                                       chr,
                                                       chromStart, chromEnd,
                                                       jcn_strand, upstr_len,
                                                       downstr_len,
                                                       chr_start + upstr_len,
                                                       chr_start + upstr_len + intron_len - 1,
                                                       intron_info[2])
                else:
                    jcn2JcnInfo[jcn_str] = JcnInfo(jcn_str,
                                                   chr,
                                                   chromStart, chromEnd,
                                                   jcn_strand, upstr_len,
                                                   downstr_len,
                                                   chr_start + upstr_len,
                                                   chr_start + upstr_len + intron_len - 1,
                                                   intron_info[2])

                if known_junctions:
                    jcn2type[jcn_str] = type

    return jcn2JcnInfo, jcn2type
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
