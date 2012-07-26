#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# combine_clusterASExons2_by_chr.py
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
    opt_parser.add_option("--root_dir",
                          dest="root_dir",
                          type="string",
                          help="""Root directory that contains each sample
                                  output files and subfolders for each
                                  chromosome.""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_file",
                          type="string",
                          help="""Output file that will contain all exclusion
                                  inclusion counts for every sample.""",
                          default=None)
    opt_parser.add_option("--left_intron",
                          dest="left_intron_file",
                          type="string",
                          help="""Output file that will contain the left side
                                  of intron retention events.  Significant
                                  intron retention events are identified in a
                                  different way.""",
                          default=None)
    opt_parser.add_option("--right_intron",
                          dest="right_intron_file",
                          type="string",
                          help="""Output file that will contain the right side
                                  of intron retention events.  Significant
                                  intron retention events are identified in a
                                  different way.""",
                          default=None)
    opt_parser.add_option("--remove_tmp_files",
                          dest="remove_tmp_files",
                          type="string",
                          help="""Remove temporary files.""",
                          default=None)
    opt_parser.add_option("--remove_interm_files",
                          dest="remove_interm_files",
                          type="string",
                          help="""Removes intermediate output files from
                                  individual getASEventReadCount runs of the sample and
                                  chromosomes. Use only to save disk space.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--root_dir")
    opt_parser.check_required("-o")
    opt_parser.check_required("--left_intron")
    opt_parser.check_required("--right_intron")

    root_dir = formatDir(options.root_dir)

    output_file = open(options.output_file,"w")
    left_intron_file = open(options.left_intron_file, "w")
    right_intron_file = open(options.right_intron_file, "w")

    # Get chr list from tmp files that are in the root dir
    chr_list = getChrs(root_dir)

    # Combine files
    for this_chr in chr_list:
        # Add main output file lines         
        
        outfile = root_dir + "/tmp_clusterASExons2_%s.out" % this_chr
        chr_out = open(outfile)
        for line in chr_out:
            output_file.write(line)
        chr_out.close()       
    
        # Add left_intron output file lines         
        left_file = root_dir + "/tmp_clusterASExons2_%s_left_intron.out " % this_chr
        chr_left_intron = open(left_file)
        for line in chr_left_intron:
            left_intron_file.write(line)
        chr_left_intron.close()

        # Add right_intron output file lines         
        right_file = root_dir + "/tmp_clusterASExons2_%s_right_intron.out " % this_chr
        chr_right_intron = open(right_file)
        for line in chr_right_intron:
            right_intron_file.write(line)
        chr_right_intron.close()

    if options.remove_tmp_files:
        os.remove(outfile)
        os.remove(left_file)
        os.remove(right_file)

    # Delete intermediate files
    if options.remove_interm_files:
        for samp in os.listdir(root_dir):
            for this_chr in chr_list:
                sub_dir = "%s/%s/%s_%s" % (root_dir,
                                           samp,
                                           samp, this_chr)   
                for filename in os.listdir(sub_dir):
                    os.remove(sub_dir + "/" + filename)
                # Write a message
                readme_file = open(sub_dir + "/README", "w")
                readme_file.write("Files deleted on %s\n" % time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime()))
                readme_file.close()
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
            if "clusterASExons2" in subfile:
                first_split = subfile.split("_")
                chr_list.append(first_split[-1].replace(".out",""))

    return chr_list

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
