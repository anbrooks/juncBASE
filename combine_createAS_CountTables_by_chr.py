#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
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
                          dest="output_file_suffix",
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
    opt_parser.add_option("--remove_tmp_files",
                          dest="remove_tmp_files",
                          action="store_true",
                          help="""Remove temporary files.""",
                          default=None)
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

    suffix = options.output_file_suffix

    output_file = open(suffix + "_AS_exclusion_inclusion_counts.txt" ,"w")
    left_intron_file = open(suffix + "_left_intron_counts.txt", "w")
    right_intron_file = open(suffix + "_right_intron_counts.txt", "w")

    lenNorm_output_file = open(suffix + "_AS_exclusion_inclusion_counts_lenNorm.txt" ,"w")
    lenNorm_left_intron_file = open(suffix + "_left_intron_counts_lenNorm.txt", "w")
    lenNorm_right_intron_file = open(suffix + "_right_intron_counts_lenNorm.txt", "w")

    # Get chr list from tmp files that are in the root dir
    chr_list = getChrs(root_dir)

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

        for line in chr_out:
            if line.startswith("#"):
                continue
            output_file.write(line)
        chr_out.close()       
       
        for line in chr_out_lenNorm:
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
            
        for line in chr_left_intron:
            if line.startswith("#"):
                continue
            left_intron_file.write(line)
        chr_left_intron.close()

        for line in chr_left_intron_lenNorm:
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
            
        for line in chr_right_intron:
            if line.startswith("#"):
                continue
            right_intron_file.write(line)
        chr_right_intron.close()

        for line in chr_right_intron_lenNorm:
            if line.startswith("#"):
                continue
            lenNorm_right_intron_file.write(line)
        chr_right_intron_lenNorm.close()

    if options.remove_tmp_files:
        os.remove(outfile)
        os.remove(left_file)
        os.remove(right_file)

        os.remove(outfile_lenNorm)
        os.remove(left_file_lenNorm)
        os.remove(right_file_lenNorm)

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
    chr_list = []

    for subfile in os.listdir(root_dir):
        if not os.path.isdir(root_dir + "/" + subfile):
            if "intron" in subfile:
                continue
            if "createAS_CountTables" in subfile:
                first_split = subfile.split("_")
                chr_list.append(first_split[3])
#                chr_list.append(first_split[-1].replace(".out",""))

    return chr_list

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
