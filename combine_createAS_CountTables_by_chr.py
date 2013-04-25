#!/usr/bin/env python
# combine_createAS_CountTables_by_chr.py
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb

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
    opt_parser.add_option("-d",
                          dest="root_dir",
                          type="string",
                          help="""Root directory of folder specified in the -d
                                  option of run_createAS_CountTables.py""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_file_prefix",
                          type="string",
                          help="""Output files that will contain all exclusion
                                  inclusion counts for every sample. As well as
                                  files for intron retention calculation.
                                  Finally, length-normalized counts are also
                                  produced.""",
                          default=None)
#   opt_parser.add_option("-o",
#                         dest="output_file",
#                         type="string",
#                         help="""Output file that will contain all exclusion
#                                 inclusion counts for every sample.""",
#                         default=None)
#   opt_parser.add_option("--left_intron",
#                         dest="left_intron_file",
#                         type="string",
#                         help="""Output file that will contain the left side
#                                 of intron retention events.  Significant
#                                 intron retention events are identified in a
#                                 different way.""",
#                         default=None)
#   opt_parser.add_option("--right_intron",
#                         dest="right_intron_file",
#                         type="string",
#                         help="""Output file that will contain the right side
#                                 of intron retention events.  Significant
#                                 intron retention events are identified in a
#                                 different way.""",
#                         default=None)
    opt_parser.add_option("--keep_tmp_files",
                          dest="keep_tmp_files",
                          action="store_true",
                          help="""Temporary files will be kept and not deleted""",
                          default=False)
#   opt_parser.add_option("--remove_interm_files",
#                         dest="remove_interm_files",
#                         type="string",
#                         help="""Removes intermediate output files from
#                                 individual getASEventReadCount runs of the sample and
#                                 chromosomes. Use only to save disk space.""",
#                         default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-d")
    opt_parser.check_required("-o")
#   opt_parser.check_required("--left_intron")
#   opt_parser.check_required("--right_intron")

    root_dir = formatDir(options.root_dir)

    prefix = options.output_file_prefix

    output_file = open(prefix + "_AS_exclusion_inclusion_counts.txt" ,"w")
    left_intron_file = open(prefix + "_left_intron_counts.txt", "w")
    right_intron_file = open(prefix + "_right_intron_counts.txt", "w")

    lenNorm_output_file = open(prefix + "_AS_exclusion_inclusion_counts_lenNorm.txt" ,"w")
    lenNorm_left_intron_file = open(prefix + "_left_intron_counts_lenNorm.txt", "w")
    lenNorm_right_intron_file = open(prefix + "_right_intron_counts_lenNorm.txt", "w")

    # Get chr list from tmp files that are in the root dir
    chr_list = list(getChrs(root_dir))

    # Getting the header from an arbitrary chr
    header = grabHeader(root_dir + "/tmp_createAS_CountTables_%s_AS_exclusion_inclusion_counts.txt" % chr_list[0])
    output_file.write(header)
    left_intron_file.write(header)
    right_intron_file.write(header)
    lenNorm_output_file.write(header)
    lenNorm_left_intron_file.write(header)
    lenNorm_right_intron_file.write(header)

    # Combine files
    for this_chr in chr_list:
        # Add main output file lines         
        
        outfile = root_dir + "/tmp_createAS_CountTables_%s_AS_exclusion_inclusion_counts.txt" % this_chr
        outfile_lenNorm = root_dir + "/tmp_createAS_CountTables_%s_AS_exclusion_inclusion_counts_lenNorm.txt" % this_chr
        try:
            chr_out = open(outfile)
            chr_out_lenNorm = open(outfile_lenNorm)
        except:
            print "No file for %s" % this_chr
            continue

        # Check if there is anything in the file
        chr_out_lines = chr_out.readlines()
        chr_out_lenNorm_lines = chr_out_lenNorm.readlines()

        if (len(chr_out.readlines()) == 1) or (len(chr_out_lenNorm_lines) == 1):
            print "Warning: No information for %s" % this_chr
            print "This may have failed in previous step."
            continue

        for line in chr_out_lines:
            if line.startswith("#"):
                continue
            output_file.write(line)
        chr_out.close()       
       
        for line in chr_out_lenNorm_lines:
            if line.startswith("#"):
                continue
            lenNorm_output_file.write(line)
        chr_out_lenNorm.close()
    
        # Add left_intron output file lines         
        left_file = root_dir + "/tmp_createAS_CountTables_%s_left_intron_counts.txt" % this_chr
        left_file_lenNorm = root_dir + "/tmp_createAS_CountTables_%s_left_intron_counts_lenNorm.txt" % this_chr
        try:
            chr_left_intron = open(left_file)
            chr_left_intron_lenNorm = open(left_file_lenNorm)
        except:
            print "No left file for %s" % this_chr
            continue

        chr_left_intron_lines = chr_left_intron.readlines()
        chr_left_intron_lenNorm_lines = chr_left_intron_lenNorm.readlines()

        if (len(chr_left_intron_lines) == 1) or (len(chr_left_intron_lenNorm_lines) == 1):
            print "Warning: No information for left intron %s" % this_chr
            print "This may have failed in previous step."
            continue
            
        for line in chr_left_intron_lines:
            if line.startswith("#"):
                continue
            left_intron_file.write(line)
        chr_left_intron.close()

        for line in chr_left_intron_lenNorm_lines:
            if line.startswith("#"):
                continue
            lenNorm_left_intron_file.write(line)
        chr_left_intron_lenNorm.close()

        # Add right_intron output file lines         
        right_file = root_dir + "/tmp_createAS_CountTables_%s_right_intron_counts.txt" % this_chr
        right_file_lenNorm = root_dir + "/tmp_createAS_CountTables_%s_right_intron_counts_lenNorm.txt" % this_chr
        try:
            chr_right_intron = open(right_file)
            chr_right_intron_lenNorm = open(right_file_lenNorm)
        except:
            print "No right file for %s" % this_chr
            continue 

        chr_right_intron_lines = chr_right_intron.readlines()
        chr_right_intron_lenNorm_lines = chr_right_intron_lenNorm.readlines()

        if (len(chr_right_intron_lines) == 1) or (len(chr_right_intron_lenNorm_lines) == 1):
            print "Warning: No information for right intron %s" % this_chr
            print "This may have failed in previous step."
            continue
            
        for line in chr_right_intron_lines:
            if line.startswith("#"):
                continue
            right_intron_file.write(line)
        chr_right_intron.close()

        for line in chr_right_intron_lenNorm_lines:
            if line.startswith("#"):
                continue
            lenNorm_right_intron_file.write(line)
        chr_right_intron_lenNorm.close()

        if not options.keep_tmp_files:
            os.remove(outfile)
            os.remove(left_file)
            os.remove(right_file)

            os.remove(outfile_lenNorm)
            os.remove(left_file_lenNorm)
            os.remove(right_file_lenNorm)

    output_file.close() 
    left_intron_file.close() 
    right_intron_file.close() 

    lenNorm_output_file.close() 
    lenNorm_left_intron_file.close()
    lenNorm_right_intron_file.close() 

#   # Delete intermediate files
#   if options.remove_interm_files:
#       for samp in os.listdir(root_dir):
#           for this_chr in chr_list:
#               sub_dir = "%s/%s/%s_%s" % (root_dir,
#                                          samp,
#                                          samp, this_chr)   
#               for filename in os.listdir(sub_dir):
#                   os.remove(sub_dir + "/" + filename)
#               # Write a message
#               readme_file = open(sub_dir + "/README", "w")
#               readme_file.write("Files deleted on %s\n" % time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime()))
#               readme_file.close()
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getChrs(root_dir):
    chr_set = set([])

    for subfile in os.listdir(root_dir):
        if not os.path.isdir(root_dir + "/" + subfile):
            if "intron" in subfile:
                continue
            if "createAS_CountTables" in subfile:
                first_split = subfile.split("_")
                chr_set.add(first_split[3])

    return chr_set

def grabHeader(file_name):
    file = open(file_name)
    header = None
    for line in file:
        if line.startswith("#"):
            header = line
            break
    file.close()
    return header
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
