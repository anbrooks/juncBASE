#!/lab/64/bin/python
# intron_exon_jctn_counts.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Counts the number of reads hit to an exon-intron boundary, above a certain offset cutoff.
"""
import os
import sys
import optparse 
import pdb
import random

from helperFunctions import runCmd,updateDictOfLists
from getASEventReadCounts import convertCoordStr
#############
# CONSTANTS #
#############
BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
COUNTING_SCRIPT = "%s/coordReadCounts.py" % BIN_DIR
if not os.path.exists(COUNTING_SCRIPT):
    print "ERROR: coordReadCounts.py needs to be in the same directory."
    sys.exit(1)

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
                          dest="intron_coords",
                          type="string",
                          help="""File of intron coordinates.  Format:
                                  type, chr, strand, start, end""",
                          default=None)
    opt_parser.add_option("-b",
                          dest="bed_intron_coords",
                          type="string",
                          help="BED file of intron coordinates.",
                          default=None)
    opt_parser.add_option("-a",
                          dest="read_alignments",
                          type="string",
                          help="""File of alignments to genome. 
                                  Format:
                                  chr, start, strand""",
                          default=None)
    opt_parser.add_option("-f",
                          dest="flanking_dist",
                          type="int",
                          help="""Distance away from exon intron junction to
                                  check for reads in.""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="offsets",
                          type="int",
                          help="""Minimum number of offsets required at each
                                  exon/intron junction. Default=1""",
                          default=1)
    opt_parser.add_option("-l",
                          dest="read_length",
                          type="int",
                          help="Length of the reads.",
                          default=1)
    opt_parser.add_option("--out_dir",
                          dest="out_dir",
                          type="string",
                          help="Output files are put here.",
                          default=None)
    opt_parser.add_option("--out_prefix",
                          dest="prefix",
                          type="string",
                          help="Prefix attached to all output files.",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-a")
    opt_parser.check_required("-f")
    opt_parser.check_required("-l")
    opt_parser.check_required("--out_dir")
    opt_parser.check_required("--out_prefix")

    # Check that the COUNTING_SCRIPT path is valid
    if not os.path.exists(COUNTING_SCRIPT):
        print "Please change COUNTING_SCRIPT path."
        opt_parser.print_help()
        sys.exit(1)

    if options.intron_coords and options.bed_intron_coords:
        print "Only one type of intron coord can be used as input." 
        opt_parser.print_help()
        sys.exit(1)

    if (not options.intron_coords) and (not options.bed_intron_coords):   
        print " Need to specify intron coordinates. See options -i or -b"
        opt_parser.print_help()
        sys.exit(1)

    intron_coords = None
    isBedFormat = False
    if options.intron_coords:
        intron_coords = open(options.intron_coords)
    if options.bed_intron_coords:
        intron_coords = open(options.bed_intron_coords)
        isBedFormat = True
        
    read_alignments = options.read_alignments

    read_length = options.read_length

    flanking_dist = options.flanking_dist
    offsets = options.offsets

    prefix = options.prefix
    out_dir = options.out_dir

    if not out_dir.endswith("/"):
        out_dir += "/"

    if not os.path.exists(out_dir):
        print "Output directory does not exist"
        sys.exit(1)

    # Intermediate Output Files
    out_coords_file = out_dir + prefix + "_intron_exon_junction_coords.out"
    out_coords = open(out_coords_file, "w")

    out_read_assoc_file = out_dir + prefix + "_intron_exon_junction_coords_w_read.out"
    
    # Final output
    out_file_name = out_dir + prefix + "_intron_exon_junction_counts.txt"
    out_file = open(out_file_name, "w")

    confident_ie_name = out_dir + prefix + "_confident_ie.txt"
    confident_ie_file = open(confident_ie_name, "w")
   
    # {intron_coord: {"left": (chr, start, end),
    #                 "right": (chr, start, end)}
    # "left" and "right" being the region at the left or right side of the
    # junction, around the exon/intron junction
    # The dict is the above but reverse mapping
    left_region_coord2intron = {}     
    right_region_coord2intron = {}     

    # {intron_coord_str:{"left":{pos:count}, 
    #                    "right":{pos:count}}
    intron_dict = {}

    regions_set = set([])

    for line in intron_coords:
        line = formatLine(line)

        if isBedFormat:
            if line.startswith("track"):
                continue
            chr, start_str, end_str = parseBEDLine(line)
        else:
            type, chr, strand, start_str, end_str = line.split("\t")

        if chr.startswith("chr"):
            chr = chr.replace("chr", "")

        intron_coord_str = "%s:%s-%s" % (chr, start_str, end_str)

        if intron_coord_str not in intron_dict:
            intron_dict[intron_coord_str] = {"left": {},
                                             "right": {}}

        start = int(start_str)
        end = int(end_str)

        left_coord = (chr,
                      start - flanking_dist,
                      start + flanking_dist - 1)

        right_coord = (chr,
                       end - flanking_dist + 1,
                       end + flanking_dist)

        updateDictOfLists(left_region_coord2intron, left_coord,
                          intron_coord_str)
        updateDictOfLists(right_region_coord2intron, right_coord,
                          intron_coord_str)

        regions_set.add(left_coord)
        regions_set.add(right_coord)

    # Print out regions out_coords
    for region_coord in regions_set:

        out_line = "%s\t%d\t%d\n" % (region_coord[0],
                                     region_coord[1],
                                     region_coord[2])

        out_coords.write(out_line)

    out_coords.close()

    # Used to make unique name for tmp file in case a shared directory is being
    # used for runs.
    rand_num = random.randrange(1,100000)

    # Get Read Counts
    print "Getting Counts in Region"
    cmd = "python %s --reads %s -l %d --coords %s -o %stmp%d.txt --read_assoc %s" % (COUNTING_SCRIPT,
                                                                          read_alignments,
                                                                          read_length,
                                                                          out_coords_file,
                                                                          out_dir,
                                                                          rand_num,
                                                                          out_read_assoc_file)
    print cmd
#    runCmd(cmd, SHELL)
    os.system(cmd)

    # Remove the tmp file
#    runCmd("rm %stmp%d.txt" % (out_dir, rand_num), SHELL)
    os.system("rm %stmp%d.txt" % (out_dir, rand_num))

    print "Getting Left and Right Counts"
    # Parse read_assoc_file to get information
    read_assoc_file = open(out_read_assoc_file)

    for line in read_assoc_file:
        line = formatLine(line)

        line_list = line.split("\t")

        read_start, read_end = getReadStartEnd(line_list[1]) 

        region_coord = getRegionCoord(line_list[2])
        intron_coord_list = getIntronStartEnds(left_region_coord2intron,
                                               right_region_coord2intron,
                                               region_coord)

       
        if region_coord in left_region_coord2intron: 
            for intron_str in left_region_coord2intron[region_coord]:
                # Put in left dictionaries
                if read_end not in intron_dict[intron_str]["left"]:
                    intron_dict[intron_str]["left"][read_end] = 1
                else:
                    intron_dict[intron_str]["left"][read_end] += 1

        if region_coord in right_region_coord2intron:
            for intron_str in right_region_coord2intron[region_coord]:
                # Check right dictionary
                if read_end not in intron_dict[intron_str]["right"]:
                    intron_dict[intron_str]["right"][read_end] = 1
                else:
                    intron_dict[intron_str]["right"][read_end] += 1

    # Print output
    confident_ie_set = set([])
    for intron_str in intron_dict:
#       chr, intron_start_str, intron_end_str = intron_str.split("_")
#       intron_start = int(intron_start_str)
#       intron_end = int(intron_end_str) 
        chr, intron_start, intron_end = convertCoordStr(intron_str)

        # Get left_counts
        if len(intron_dict[intron_str]["left"]) >= offsets:
            left_count = getTotalCounts(intron_dict[intron_str]["left"])
            confident_ie = "%s:%d-%d" % (chr, intron_start - 1, intron_start)
            confident_ie_set.add(confident_ie)
        else:
            left_count = 0

        # Get right counts
        if len(intron_dict[intron_str]["right"]) >= offsets:
            right_count = getTotalCounts(intron_dict[intron_str]["right"])
            confident_ie = "%s:%d-%d" % (chr, intron_end, intron_end + 1)
            confident_ie_set.add(confident_ie)
        else:
            right_count = 0

        if left_count == 0 and right_count == 0:
            continue

        print_line = "%s\t%d\t%d\n" % (intron_str,
                                     left_count,
                                     right_count)

        out_file.write(print_line)

    # Now print out confident set of ie
    for ie in confident_ie_set:
        confident_ie_file.write("%s\n" % ie)

    confident_ie_file.close()
	
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

def getIntronStartEnds(left_region_coord2intron, right_region_coord2intron,
                       region_coord):

    returnList = []

    if region_coord in left_region_coord2intron:
        returnList.extend(left_region_coord2intron[region_coord])
    if region_coord in right_region_coord2intron:
        returnList.extend(right_region_coord2intron[region_coord])    

    return returnList

def getReadStartEnd(read_str):
    """
    Returns the start and end position
    """
#    elems = read_str.split("_")
    elems = convertCoordStr(read_str)

    return elems[1], elems[2]

def getRegionCoord(region_coord_str):
    region_coord = convertCoordStr(region_coord_str)
#   region_coord_str_elems = region_coord_str.split("_")
#   region_coord = (region_coord_str_elems[0],
#                   int(region_coord_str_elems[1]),    
#                   int(region_coord_str_elems[2]))

    return region_coord

def getTotalCounts(count_dict):
    total = 0
    
    for overhang_pos in count_dict:
        total += count_dict[overhang_pos]

    return total

def parseBEDLine(line):
    """
    Returns the chr, start, and end of intron in string format.
    """
    input_list = line.split("\t")

    chr = input_list[0]

    chromStart = int(input_list[1])

    blockSizes = map(int, input_list[-2].split(","))
    blockStarts = map(int, input_list[-1].split(","))
   
    intron_start = chromStart + blockSizes[0] + 1
    intron_end = chromStart + blockStarts[1]

    return chr, repr(intron_start), repr(intron_end)
    
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
