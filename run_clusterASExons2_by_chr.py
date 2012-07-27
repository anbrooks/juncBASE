#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# run_clusterASExons2_by_chr.py
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
from broad_helperFunctions import runLSF

from createPseudoSample import getChr
#############
# CONSTANTS #
#############
DEF_NUM_PROCESSES = 2

BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/clusterASExons2.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: clusterASExons2.py needs to be in the same directory."
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
    opt_parser.add_option("-d",
                          dest="root_dir",
                          type="string",
                          help="""Root directory that contains subdirectoires
                                  with output from getASEventReadCounts""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="input_dir",
                          type="string",
                          help="""Directory containing original input files to
                                  getASEventReadCounts.py. This is used to
                                  obtain the chromosome information.""",
                          default=None)
    opt_parser.add_option("-s",
                          dest="samples",
                          type="string",
                          help="""Comma separated list of the samples that will
                                  be used.  The order which they are given is
                                  the order in the output of the file.""",
                          default=None)
    opt_parser.add_option("--lengthNorm",
                          dest="lengthNorm",
                          action="store_true",
                          help="""Flag to indicate length normalization was
                                  done on the counts. Used for splitting the IR
                                  counts back into left and right counts""",
                          default=False)
    opt_parser.add_option("--num_processes",
                          dest="num_processes",
                          type="int",
                          help="""Will run each chromosome in batches using this
                                  number of parallel processes. DEF=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)
    opt_parser.add_option("--run_LSF",
                          dest="run_lsf",
                          action="store_true",
                          help="Will run everything through LSF",
                          default=False)


    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-d")
    opt_parser.check_required("-i")
    opt_parser.check_required("-s")

    root_dir = formatDir(options.root_dir)
    # Change to the root directory to make sure output files are put here
    os.chdir(root_dir)

    input_dir = formatDir(options.input_dir)

    samples = options.samples

    lengthNorm = options.lengthNorm

    num_processes = options.num_processes
    run_lsf = options.run_lsf

    chr_list = getChr(input_dir)

    ctr = 0
    for this_chr in chr_list:
        ctr += 1

        cmd = "python %s " % SCRIPT
        cmd += "-d %s " % root_dir
        cmd += "-o tmp_clusterASExons2_%s.out " % this_chr
        cmd += "--left_intron tmp_clusterASExons2_%s_left_intron.out " % this_chr
        cmd += "--right_intron tmp_clusterASExons2_%s_right_intron.out " % this_chr
        cmd += "-s %s " % samples

        if lengthNorm:
            cmd += "--lengthNorm "

        cmd += "--which_chr %s" % this_chr

        if run_lsf:
            runLSF(cmd,
                   "%s.clusterASExons2.bsub.out" % this_chr,
                   samples.replace(",","-") + "_" + this_chr,
                   "hour")
            continue

        if ctr % num_processes == 0:
            os.system(cmd)
        else:
            print cmd
            Popen(cmd, shell=True, executable=SHELL)

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
