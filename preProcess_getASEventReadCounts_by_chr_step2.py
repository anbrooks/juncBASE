#!/lab/64/bin/python
# preProcess_getASEventReadCounts_step2.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Creates the intron_exon junction files.
"""

import sys
import optparse 
import os
import pdb

from subprocess import Popen
from helperFunctions import runCmd
from merge_by_chr_juncBASE_junction_file import isBadChr
from getASEventReadCounts import convertCoordStr, formatCoordStr
from preProcess_getASEventReadCounts import DEF_BEG_JCN_OVERHANG as DEF_JCN_OVERHANG
#############
# CONSTANTS #
#############
#DEF_NUM_PROCESSES = 1
PSEUDO_COUNT = 3
DEF_CONFIDENCE = 3

SHELL = "/bin/tcsh"
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
    opt_parser.add_option("-i",
                          dest="input_dir",
                          type="string",
                          help="""Root directory containing all input
                                  directories that will be used for
                                  getASEventReadCounts.py""",
                          default=None)
    opt_parser.add_option("--by_chr",
                          dest="by_chr",
                          action="store_true",
                          help="""Indicates that each sample has the input files
                                  split by chr. Expects the directory structure
                                  created through the first preProcess step.""",
                          default=False)
    opt_parser.add_option("--no_clean_chr",
                          action="store_true",
                          help="""By default, will only incorporate standard
                                  chromosomes. This option will allow any
                                  additional non-standard chr through.""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")

    input_dir = options.input_dir

    if os.path.exists(input_dir):
        input_dir = os.path.abspath(options.input_dir)
    else:
        print "ERROR. Input directory does not exist: %s" % input_dir
        opt_parser.print_help()
        sys.exit(1)

    if not input_dir.endswith("/"):
        input_dir += "/"

    by_chr = options.by_chr
    no_clean_chr = options.no_clean_chr

    # If by_chr, this is actually chr2intron_jcn2strand    
    intron_jcn2strand = getAllJcns(input_dir, by_chr, no_clean_chr)

    # Junciton bed will contain all junctions.
    createTmpJcnBed(input_dir, intron_jcn2strand, by_chr)

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def createTmpJcnBed(input_dir, intron_junction2strand, by_chr):
    """
    if by_chr, intron_junction2strand is really chr2intron_junction2strand
    """
    chr2jcn_bed = {}

    if by_chr:
        for this_chr in intron_junction2strand:
            chr2jcn_bed[this_chr] = "%stmp_%s_preProcess_getASEventReadCounts_step2.bed" % (input_dir, 
                                                                                            this_chr)
    else:
        chr2jcn_bed["all"] = input_dir + "tmp_preProcess_getASEventReadCounts_step2.bed"
        
    for this_chr in chr2jcn_bed:
        # Create temporary junction bed
        output_bed = open(chr2jcn_bed[this_chr], "w")

        if by_chr:
            this_intron_junction2strand = intron_junction2strand[this_chr]
        else:
            this_intron_junction2strand = intron_junction2strand
                

        # Write header
        header = "track name=\"all_jcns\"\n"
        output_bed.write(header)

        for jcn_str in this_intron_junction2strand:
            (chr, chrStart, chrEnd,
             blockLens, secondBlockStart) = parse_jcn_str(jcn_str)

            bed_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,0\t2\t%d,%d\t0,%d\n" % (chr,
                                                                                    chrStart,
                                                                                    chrEnd,
                                                                                    jcn_str,
                                                                                    PSEUDO_COUNT,
                                                                                    this_intron_junction2strand[jcn_str],
                                                                                    chrStart,
                                                                                    chrEnd,
                                                                                    blockLens,
                                                                                    blockLens,
                                                                                    secondBlockStart)
            output_bed.write(bed_line)

        output_bed.close()

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getAllJcns(input_dir, by_chr, no_clean_chr):
    # GET INITIAL EXON-EXON JUNCTION
    if by_chr:
        chr2intron_junction2strand = {}
    else:
        intron_junction2strand = {}

    for samp in os.listdir(input_dir):
        if samp.endswith(".bed") or samp.endswith(".txt") or\
           samp.endswith(".py") or samp.startswith("."):
            continue
        if not os.path.isdir(input_dir + samp):
            continue
        
        if by_chr:
            for chr_dir in os.listdir(input_dir + samp):
                if not os.path.isdir(input_dir + samp + "/" + chr_dir):
                    continue
                this_chr = getChr(chr_dir, samp)
                if not no_clean_chr:
                    if isBadChr(this_chr):
                        continue

                if this_chr not in chr2intron_junction2strand:
                    chr2intron_junction2strand[this_chr] = {}

                # Now process bed file
                bed_file = open(input_dir + samp + "/" + chr_dir + "/" + chr_dir + "_junctions.bed")
                jcn_str2strand = getJunctionStrs(bed_file)

                for jcn_str in jcn_str2strand:
                    if jcn_str in chr2intron_junction2strand[this_chr]:
                        if jcn_str2strand[jcn_str] != chr2intron_junction2strand[this_chr][jcn_str]:
                            chr2intron_junction2strand[this_chr][jcn_str] = "."
                    else:
                        chr2intron_junction2strand[this_chr][jcn_str] = jcn_str2strand[jcn_str]

                bed_file.close()
 
        else:
            # Just process one bed file
            bed_file = open(input_dir + samp + "/" + samp + "_junctions.bed")

            jcn_str2strand = getJunctionStrs(bed_file)

            for jcn_str in jcn_str2strand:
                if jcn_str in intron_junction2strand:
                    if jcn_str2strand[jcn_str] != intron_junction2strand[jcn_str]:
                        intron_junction2strand[jcn_str] = "."
                else:
                    intron_junction2strand[jcn_str] = jcn_str2strand[jcn_str]

            bed_file.close()

    if by_chr:
        return chr2intron_junction2strand

    return intron_junction2strand

def getChr(chr_dir, samp):
    chrName = chr_dir.split(samp)[-1]
    chrName = chrName.strip("_")
    return chrName

def getJunctionStrs(bed_file):

    jcn2strand = {}

    for line in bed_file:
        line = formatLine(line)

        if line.startswith("track"):
            continue

        line_list = line.split("\t")

        chr = line_list[0]
        chrStart = int(line_list[1])
    
        strand = line_list[5]
        
        blocks = map(int,line_list[10].split(","))
        block_starts = map(int,line_list[11].split(","))

        intron_start = chrStart + blocks[0] + 1
        intron_end = chrStart + block_starts[1]

        jcn_str = formatCoordStr(chr, intron_start, intron_end)

        jcn2strand[jcn_str] = strand

    return jcn2strand



def parse_jcn_str(jcn_str):
    """
        (chr, chrStart, chrEnd, strand, blockLens,
         secondBlockStart) = 
    """
    (chr, intron_start, intron_end) = convertCoordStr(jcn_str)
#   (chr, intron_start_str, intron_end_str) = jcn_str.split()

#   intron_start = int(intron_start_str)
#   intron_end = int(intron_end_str)

    chrStart = intron_start - DEF_JCN_OVERHANG - 1
    chrEnd = intron_end + DEF_JCN_OVERHANG

    blockLens = DEF_JCN_OVERHANG

    length = chrEnd - chrStart
    secondBlockStart = length - blockLens

    return chr, chrStart, chrEnd, blockLens, secondBlockStart
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
