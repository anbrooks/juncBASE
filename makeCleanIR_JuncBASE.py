#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
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
    opt_parser.add_option("--in_prefix",
                          dest="prefix",
                          type="string",
                          help="Prefix of files created to make the AS tables.",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--in_prefix")

    prefix = options.prefix

    suffixes = ["_AS_exclusion_inclusion_counts.txt",
                "_left_intron_counts.txt",
                "_right_intron_counts.txt",
                "_AS_exclusion_inclusion_counts_lenNorm.txt",
                "_left_intron_counts_lenNorm.txt",
                "_right_intron_counts_lenNorm.txt"]

    for suffix in suffixes:
        cmd = "awk -F \"\t\" \'{if (!(($1 == \"N\") && ($2 == \"intron_retention\"))) print $0}\' "
        cmd += prefix + suffix
        cmd += "> %s_IR_cleaned%s" % (prefix, suffix)
        os.system(cmd)
    

    print "New files use the following prefix:"
    print prefix + "_IR_cleaned"
			
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

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
