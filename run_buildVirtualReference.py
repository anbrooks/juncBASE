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

from subprocess import Popen
#############
# CONSTANTS #
#############
BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/buildVirtualReference_median.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: buildVirtualReference_median.py needs to be in the same directory."
    sys.exit(1)

DEF_THRESH = 25
DEF_NUM_PROCESSES = 2

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
    opt_parser.add_option("--in_prefix",
                          dest="in_prefix",
                          type="string",
                          help="""Prefix of output files created from
                                  createAS_CountTables. In createAS_CountTables,
                                  this is the -o option""",
                          default=None)
    opt_parser.add_option("--total_thresh",
                          dest="total_thresh",
                          type="int",
                          help="Threshold for total number of counts. DEF=%d" % DEF_THRESH,
                          default=DEF_THRESH)
    opt_parser.add_option("--remove_from_median",
                          dest="remove_from_median",
                          type="string",
                          help="""Comma separated list of samples that will not
                                  be used when calculating the median.""",
                          default=None)
    opt_parser.add_option("--weights",
                          dest="weights",
                          type="string",
                          help="""Comma separated list of weights given in the
                                  order of the samples in the table. Weights are
                                  used to create a weighted median. Default is
                                  equal weight for all samples.""",
                          default=None)
    opt_parser.add_option("-p",
                          dest="num_processes",
                          type="int",
                          help="""Number of processes to run simultaneously.
                                  Default=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)



    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--in_prefix")

    in_prefix = options.in_prefix
    
    total_thresh = options.total_thresh
    remove_from_median = options.remove_from_median
    weights = options.weights 

    num_processes = options.num_processes

    out_prefix = in_prefix + "_w_VR"

    file_suffix_list = ["_AS_exclusion_inclusion_counts.txt",
                   "_left_intron_counts.txt",
                   "_right_intron_counts.txt",
                   "_AS_exclusion_inclusion_counts_lenNorm.txt",
                   "_left_intron_counts_lenNorm.txt",
                   "_right_intron_counts_lenNorm.txt",]

    ctr = 0
    for file_suffix in file_suffix_list:
        cmd = "python %s " % SCRIPT
        cmd += "--input %s%s " % (in_prefix, file_suffix)
        cmd += "--output_median_ratio %s%s " % (out_prefix, file_suffix)
   
        if remove_from_median: 
            cmd += "--remove_from_median %s " % remove_from_median
        if weights:
            cmd += "--weights %s " % weights

        cmd += "--total_thresh %d" % total_thresh
        
        ctr += 1

        if ctr % num_processes == 0:
            os.system(cmd)
            print cmd
        else:
            print cmd
            Popen(cmd, shell=True, executable=SHELL) 

    # Print out the prefix that should be used as a reminder
    print "Output files were created with the following prefix: %s" % out_prefix
			
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
