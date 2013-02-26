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
from helperFunctions import runLSF

from createPseudoSample import getChr
#############
# CONSTANTS #
#############
DEF_NUM_PROCESSES = 2

BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/createAS_CountTables.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: createAS_CountTables.py needs to be in the same directory."
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
    opt_parser.add_option("--jcn_seq_len",
                          dest="jcn_seq_len",
                          type="int",
                          help="""Value used in getASEventReadCounts""", 
                          default=None)
    opt_parser.add_option("-s",
                          dest="samples",
                          type="string",
                          help="""Comma separated list of the samples that will
                                  be used or a file of sample names, one per
                                  line. The order which they are given is
                                  the order in the output of the file.""",
                          default=None)
#   opt_parser.add_option("--lengthNorm",
#                         dest="lengthNorm",
#                         action="store_true",
#                         help="""Flag to indicate length normalization was
#                                 done on the counts. Used for splitting the IR
#                                 counts back into left and right counts""",
#                         default=False)
    opt_parser.add_option("--num_processes",
                          dest="num_processes",
                          type="int",
                          help="""Will run each chromosome in batches using this
                                  number of parallel processes. DEF=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)
    opt_parser.add_option("--print_cmd",
                          dest="print_cmd",
                          action="store_true",
                          help="""Will not run any processes, but print the
                                  commands""",
                          default=False)
    opt_parser.add_option("--check",
                          dest="check",
                          action="store_true",
                          help="Will check which samples are not finished.",
                          default=False)
    opt_parser.add_option("--force",
                          dest="force",
                          action="store_true",
                          help="""By default, will only run jobs that need to be
                                  completed. This will force to run all
                                  jobs.""",
                          default=False)
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
    opt_parser.check_required("--jcn_seq_len")

    root_dir = formatDir(options.root_dir)
    input_dir = formatDir(options.input_dir)
    # Change to the root directory to make sure output files are put here
    os.chdir(root_dir)


    samples = options.samples

    jcn_seq_len = options.jcn_seq_len

#    lengthNorm = options.lengthNorm

    print_cmd = options.print_cmd
    check = options.check
    force = options.force

    num_processes = options.num_processes
    run_lsf = options.run_lsf

    chr_list = getChr(input_dir)

    ctr = 0
    for this_chr in chr_list:
        files_are_present = False
        expected_out_files = ["%s/tmp_createAS_CountTables_%s_AS_exclusion_inclusion_counts.txt" % (root_dir, this_chr),
                              "%s/tmp_createAS_CountTables_%s_left_intron_counts.txt" % (root_dir, this_chr),
                              "%s/tmp_createAS_CountTables_%s_right_intron_counts.txt" % (root_dir, this_chr),
                              "%s/tmp_createAS_CountTables_%s_AS_exclusion_inclusion_counts_lenNorm.txt" % (root_dir, this_chr),
                              "%s/tmp_createAS_CountTables_%s_left_intron_counts_lenNorm.txt" % (root_dir, this_chr),
                              "%s/tmp_createAS_CountTables_%s_right_intron_counts_lenNorm.txt" % (root_dir, this_chr)]
        try: 
            for expect_file in expected_out_files:
                if os.path.getsize(expect_file) == 0:
                    files_are_present = False
                    if check:
                        print "Cannot find files for: %s" % this_chr
                    break
                else:
                    files_are_present = True
        except:
            if check:
                print "Cannot find files for: %s" % this_chr

        if check:
            continue

        if not force:
            if files_are_present:
                continue

        ctr += 1

        cmd = "python %s " % SCRIPT
        cmd += "-d %s " % root_dir
        cmd += "-o %s/tmp_createAS_CountTables_%s " % (root_dir, this_chr)
#       cmd += "--left_intron tmp_createAS_CountTables_%s_left_intron.out " % this_chr
#       cmd += "--right_intron tmp_createAS_CountTables_%s_right_intron.out " % this_chr
        cmd += "-s %s " % samples
        cmd += "--jcn_seq_len %d " % jcn_seq_len

#       if lengthNorm:
#           cmd += "--lengthNorm "

        cmd += "--which_chr %s" % this_chr

        if print_cmd:
            print cmd
            continue

        if run_lsf:
            runLSF(cmd,
                   "%s.createAS_CountTables.bsub.out" % this_chr,
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
