#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# run_preProcess_by_chr_step1.py
# Author: Angela Brooks
# Program Completion Date:
# Description: Will run batched runs of
# preProcess_getASEventReadCounts_by_chr.py
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb
import pysam

from subprocess import Popen

from helperFunctions import runCmd
from preProcess_getASEventReadCounts_by_chr import getReferences, formatChr
#############
# CONSTANTS #
#############
SCRIPT = "~/bin/preProcess_getASEventReadCounts_by_chr.py"

DEF_NUM_PROCESSES = 1
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
                          dest="input_file",
                          type="string",
                          help="""File containing sample name with bam file
                                  location.""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_dir",
                          type="string",
                          help="Root directory for output files.",
                          default=None)
    opt_parser.add_option("--preProcess_options",
                          dest="preProcess_options",
                          type="string",
                          help="""Options to pass to
                                  preProcess_getASEventReadCounts_by_chr.py.
                                  Given in quotes (\")""",
                          default=None)
    opt_parser.add_option("-p",
                          dest="num_processes",
                          type="int",
                          help="""Will batch the runs with this number of
                                  processes. DEF=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)
    opt_parser.add_option("--nice",
                          dest="nice",
                          action="store_true",
                          help="Will run processes with nice.",
                          default=False)
    opt_parser.add_option("--check",
                          dest="check",
                          action="store_true",
                          help="""Will check samples that are not done and print
                                  out which need to still be run""",
                         default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")

    input_file = open(options.input_file)

    output_dir = options.output_dir

    preProcess_options = options.preProcess_options
    num_processes = options.num_processes

    check = options.check

    ctr = 0
    for line in input_file:
        line = formatLine(line)

        samp, bam = line.split("\t")
    
        if check:
            bam = pysam.Samfile(bam, "rb")
            chr_names_unformatted = getReferences([bam])
            bam.close()

            chr_names = []
            for this_chr in chr_names_unformatted:
                chr_names.append(formatChr(this_chr))

            for chr in chr_names:
                check_file = "%s/%s/%s_%s/%s_%s_junctions.bed" % (output_dir,
                                                                  samp,
                                                                  samp, chr,
                                                                  samp, chr)
                if not os.path.exists(check_file):
                    print "Cannot find: %s" % check_file 
                    continue

                if os.path.getsize(check_file) == 0:
                    print "File is empty: %s" % check_file

            continue

        cmd = ""
        if options.nice:
            cmd = "nice "
        cmd += "python %s " % SCRIPT
        cmd += "-s %s " % bam
        cmd += "-n %s " % samp
        cmd += "-o %s " % output_dir
        cmd += "%s " % preProcess_options.strip()
        cmd += "--by_chr" 

        ctr += 1

        if ctr % num_processes == 0:
            print cmd
            runCmd(cmd, SHELL, True)
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
