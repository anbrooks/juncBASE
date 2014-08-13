#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# merge_by_chr_juncBASE_junction_file.py 
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
    opt_parser.add_option("--juncBASE_input_root",
                          dest="root_dir",
                          type="string",
                          help="""Root directory containing subdirectory of
                                  juncBASE input files. This should be the same
                                  directory used in the -o option of
                                  preProcess_getASEventReadCounts_by_chr.py""",
                          default=None)
    opt_parser.add_option("--clean_chr",
                          dest="clean_chr",
                          action="store_true",
                          help="""Will only combine standard chromosomes. e.g.,
                                  it will ignore chromosomes with GL, NC, M in
                                  the name.""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--juncBASE_input_root")

    root_dir = formatDir(options.root_dir)

    for subdir in os.listdir(root_dir):
        samp_dir = root_dir + "/" + subdir
        if not os.path.isdir(samp_dir):
            continue

        if "pseudo" in samp_dir:
            continue

        # Create junction bed file
        combined_bed = open("%s/%s_junctions.bed" % (samp_dir,
                                                     subdir),
                            "w")
        # write header
        bed_header = "track name=\"%s_jcn_counts\" itemRgb=\"On\" useScore=1\n" % subdir
        combined_bed.write(bed_header)

        for sub_subdir in os.listdir(samp_dir):

            if options.clean_chr:
                if isBadChr(sub_subdir):
                    continue

            chr_dir = samp_dir + "/" + sub_subdir
            if "pseudo" in chr_dir:
                continue

            if not os.path.isdir(chr_dir):
                continue

            bed_file = open("%s/%s_junctions.bed" % (chr_dir,
                                                     sub_subdir))
            bed_lines = bed_file.readlines()
            bed_file.close()

            bed_lines.pop(0) # remove header
            combined_bed.writelines(bed_lines)

        combined_bed.close()
        
			
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

def isBadChr(subdir):
    if "GL" in subdir:
        return True
    if "M" in subdir:
        return True
    if "NC" in subdir:
        return True
    if "random" in subdir:
        return True

    return False
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
