#!/lab/64/bin/python
# preProcess_getASEventReadCounts_step3.py
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
import gzip

from subprocess import Popen
from helperFunctions import runCmd, updateDictOfLists
#############
# CONSTANTS #
#############
#DEF_NUM_PROCESSES = 1
PSEUDO_COUNT = 3
DEF_CONFIDENCE = 3

BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
IE_SCRIPT = "%s/intron_exon_jctn_counts.py" % BIN_DIR
if not os.path.exists(IE_SCRIPT):
    print "ERROR: intron_exon_jctn_counts.py needs to be in the same directory."
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
                          dest="input_dir",
                          type="string",
                          help="""Root directory containing all input
                                  directories that will be used for
                                  getASEventReadCounts.py""",
                          default=None)
    opt_parser.add_option("-n",
                          dest="name",
                          type="string",
                          help="""Sample name.  Must be the same name as a
                                  subdirectory in the input_dir.""",
                          default=None)
    opt_parser.add_option("-t",
                          dest="tmp_jcn_file",
                          type="string",
                          help="""Location of the temporary junction file
                                  created from step 2""",
                          default=None)
    opt_parser.add_option("-l",
                          dest="read_length",
                          type="int",
                          help="""Expected read length if all reads should be of
                                  the same length""",
                          default=None)
    opt_parser.add_option("--min_overhang",
                          dest="min_overhang",
                          type="int",
                          help="""Minimum overhang used to determine a junction
                                  alignment (Used as input to find intron-exon
                                  junctions.""",
                          default=None)
    opt_parser.add_option("-c",
                          dest="confidence_score",
                          type="int",
                          help="""The mininmum number of offsets a junction
                                  has to have in order to be considered
                                  confident. Default=%s""" % DEF_CONFIDENCE,
                          default=DEF_CONFIDENCE)
#   opt_parser.add_option("--single",
#                         dest="single_reads",
#                         action="store_true",
#                         help="""Flag to indicate that there are single reads.
#                                 If no flag is given, it defaults to single
#                                 reads.""",
#                         default=False)
#   opt_parser.add_option("--paired",
#                         dest="paired_reads",
#                         action="store_true",
#                         help="""Flag to indicate that there are paired-end
#                                 reads.""",
#                         default=False)
#   opt_parser.add_option("--mixed",
#                         dest="mixed_reads",
#                         action="store_true",
#                         help="""Flag to indicate that there are both paired-end
#                                 and single reads.""",
#                         default=False)


    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-n")
    opt_parser.check_required("-t")
    opt_parser.check_required("--min_overhang")

    input_dir = options.input_dir

    forced_read_len = options.read_length
    min_overhang = options.min_overhang    

    confidence_score = options.confidence_score 
    samp = options.name

    single_reads = True
    paired_reads = False
    mixed_reads = False
#   single_reads = options.single_reads
#   paired_reads = options.paired_reads
#   mixed_reads = options.mixed_reads

#   if mixed_reads:
#       single_reads = True
#       paired_reads = True
#   elif ((not mixed_reads) and (not paired_reads) and
#         (not single_reads)):
#       single_reads = True


    tmp_jcn_file = options.tmp_jcn_file

    if os.path.exists(input_dir):
        input_dir = os.path.abspath(options.input_dir)
    else:
        print "Error. Input directory does not exists: %s" % input_dir
        opt_parser.print_help()
        sys.exit(1)

    if not input_dir.endswith("/"):
        input_dir += "/"

    if not os.path.exists(input_dir + samp):
        print "Expected path does not exist: %s" % (input_dir + samp)
        opt_parser.print_help()
        sys.exit(1)

    # Now run IE script with this junction file and all other sample genome
    # files
    if single_reads:
        genome_file_name = input_dir + samp + "/" + samp + "_genome_reads.txt.gz"
    if paired_reads:
        paired_genome_file_name = input_dir + samp + "/" + samp + "_paired_end_as_single_genome.txt"

    len2readFile = {}
    readFile2handle = {}

    if paired_reads:
        if os.path.exists(paired_genome_file_name):
            paired_read_set = update_len2readFile(len2readFile, readFile2handle, 
                                                  input_dir, samp, paired_genome_file_name, True)
        else:
            print "Cannot find paired end file: %s" % paired_genome_file_name
            sys.exit(1)

    if single_reads:
        update_len2readFile(len2readFile, readFile2handle, input_dir, samp, genome_file_name) 

    # Close all filehandles
    for read_file in readFile2handle:
        readFile2handle[read_file].close()

    if forced_read_len:
        if len(len2readFile) > 1:
            print "More than one read length in reads."
            opt_parser.print_help()    
            sys.exit(1)
        if forced_read_len != len2readFile.keys()[0]:
            print "Read length in file (%d) does not match specified read length (%d)." % (len2readFile.keys()[0],
                                                                                           forced_read_len)
            opt_parser.print_help()    
            sys.exit(1)

    for read_len in len2readFile:
        # Now run the script to get the intron/exon junction counts
        cmd = "python %s " % IE_SCRIPT
        cmd += "-b %s " % tmp_jcn_file
        cmd += "-a %s " % len2readFile[read_len]
        cmd += "--out_dir %s " % (input_dir + samp + "/")
        cmd += "--out_prefix %s " % ("tmp_" + samp + "_" + repr(read_len))
        cmd += "-l %d " % read_len            
        cmd += "-f %d " % (read_len - min_overhang)
        cmd += "-o %d" % confidence_score

#        runCmd(cmd, SHELL, True)
        os.system(cmd)

    # Combine all count files into one sum (*_intron_exon_junction_counts.txt)
    all_introns = combineCounts(input_dir, len2readFile.keys(), samp)
    
    # Find paired reads and create file associating paired reads.
    if paired_reads:
        make_paired_end_ie_junctions2qname(all_introns, input_dir, len2readFile.keys(),
                                           min_overhang, samp, paired_read_set)


		
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def combineCounts(input_dir, read_lengths, samp):
    new_count_file_name = input_dir + samp + "/" + samp + "_intron_exon_junction_counts.txt"
    new_count_file = open(new_count_file_name, "w")
    
    intron2counts = {}

    for read_len in read_lengths:
        len_count_file_name = input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_intron_exon_junction_counts.txt"
        len_count_file = open(len_count_file_name)

        for line in len_count_file:
            line = formatLine(line)

            line_list = line.split("\t")

            intron = line_list[0]
            left_count = int(line_list[1])
            right_count = int(line_list[2])

            if intron in intron2counts:
                intron2counts[intron][0] += left_count
                intron2counts[intron][1] += right_count
            else:
                intron2counts[intron] = [left_count, right_count] 

        len_count_file.close()

        # Remove temporary files
        rm_cmd = "rm " + input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_intron_exon_junction_counts.txt"
        os.system(rm_cmd)

        rm_cmd = "rm " + input_dir + samp + "/tmp_" + repr(read_len) + "_reads.txt.gz"
#        runCmd(rm_cmd, SHELL, True)
        os.system(rm_cmd)

        rm_cmd = "rm " + input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_intron_exon_junction_coords_w_read.out"
#        runCmd(rm_cmd, SHELL, True)
        os.system(rm_cmd)

        rm_cmd = "rm " + input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_intron_exon_junction_coords.out"
#        runCmd(rm_cmd, SHELL, True)
        os.system(rm_cmd)

        rm_cmd = "rm " + input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_confident_ie.txt"
        os.system(rm_cmd)

    # Now print results
    for intron in intron2counts:  
        outline = "%s\t%d\t%d\n" % (intron,
                                    intron2counts[intron][0],
                                    intron2counts[intron][1])    

        new_count_file.write(outline)

    new_count_file.close()

    return intron2counts.keys()

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line


def getTheseRegionCoords(all_introns, read_length, overhang):
    """
    Given the read length and the overhang, calculate the region coord around
    the splice sites of the intron
    """
    flanking_dist = read_length - overhang

    region_coord2ie = {}

    for intron in all_introns:
        chr, start_str, end_str = intron.split("_")

        intron_start = int(start_str)
        intron_end = int(end_str)

        left_coord = "%s_%d_%d" % (chr,
                                   intron_start - flanking_dist,
                                   intron_start + flanking_dist - 1)
        right_coord = "%s_%d_%d" % (chr,
                                    intron_end - flanking_dist + 1,
                                    intron_end + flanking_dist)        

        region_coord2ie[left_coord] = "%s_%d_%d" % (chr,
                                                    intron_start - 1,
                                                    intron_start)
        region_coord2ie[right_coord] = "%s_%d_%d" % (chr,
                                                     intron_end,
                                                     intron_end + 1)

    return region_coord2ie
           

def make_paired_end_ie_junctions2qname(all_introns, input_dir, read_lengths, overhang, samp, paired_read_set):

    ie2qname_file_name = input_dir + samp + "/" + samp + "_paired_end_ie_junctions2qname.txt"
    ie2qname_file = open(ie2qname_file_name, "w")

    ie2qnames = {}

    for read_len in read_lengths:
        # Set of region coords that will be searched within the region and read
        # association file
        region_coord2ie = getTheseRegionCoords(all_introns, read_len,
                                                overhang) 

        confident_ie_file_name = input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_confident_ie.txt"
        confident_ie_file = open(confident_ie_file_name)

        confident_ie_set = set([])
        for line in confident_ie_file:
            confident_ie_set.add(formatLine(line))
        confident_ie_file.close()

        coords_w_reads_file_name = input_dir + samp + "/tmp_" + samp + "_" + repr(read_len) + "_intron_exon_junction_coords_w_read.out"
        coords_w_reads_file = open(coords_w_reads_file_name)
        
        for line in coords_w_reads_file:
            line = formatLine(line)
            line_list = line.split("\t")

            if line_list[-1] in region_coord2ie:
                ie = region_coord2ie[line_list[-1]]
                if ie in confident_ie_set:
                    updateDictOfLists(ie2qnames, ie, line_list[0])                

        coords_w_reads_file.close()           
   
    for ie in ie2qnames:
        outline = "%s\t%s\n" % (ie,
                              ",".join(ie2qnames[ie]))     
        ie2qname_file.write(outline)

    ie2qname_file.close() 

def update_len2readFile(len2readFile, readFile2handle, input_dir, samp, genome_file_name, isPaired=False):
    """
    Returns set of paired_end reads qnames 
    """
    paired_read_set = None
    if isPaired:
        paired_read_set = set([])

    genome_file = gzip.open(genome_file_name, "rb")

    for line in genome_file:
        line = formatLine(line)
        
        line_list = line.split("\t")

        read_length = int(line_list[-1]) - int(line_list[-2]) + 1       

        if read_length in len2readFile:
            readFile2handle[len2readFile[read_length]].write(line + "\n")
        else:
            readFile_name = input_dir + samp + "/tmp_" + repr(read_length) + "_reads.txt.gz"
            len2readFile[read_length] = readFile_name
            readFile2handle[readFile_name] = gzip.open(readFile_name, "wb")

            readFile2handle[readFile_name].write(line + "\n")
        if isPaired:
            paired_read_set.add(line_list[0])



    genome_file.close()

    return paired_read_set
    
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
