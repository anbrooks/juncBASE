#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# runCufflinks.py
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os

from helperFunctions import runCmd
from subprocess import Popen
from broad_helperFunctions import runLSF
#############
# CONSTANTS #
#############
SHELL = "/bin/tcsh"

CUFF_EXEC = "/home/unix/brooks/src/Cufflinks/cur/cufflinks"

DEF_NUM_PROCESSES = 2
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
                          dest="input",
                          type="string",
                          help="""Tab-delimited file that specifies sample name
                                  and bam location""",
                          default=None)
    opt_parser.add_option("--force",
                          dest="force",
                          action="store_true",
                          help="""By default, will only run Cufflinks if no
                                  output file exists. This option forces the
                                  runs on every sample.""",
                          default=False)
    opt_parser.add_option("--txt_ref",
                          dest="txt_ref",
                          type="string",
                          help="Transcript reference used for assembly.",
                          default=None)
    opt_parser.add_option("--quantitate",
                          dest="quantitate",
                          action="store_true",
                          help="""Will quantitate against reference transcript
                                  annotations instead of assembly""",
                          default=False)
    opt_parser.add_option("--out_dir",
                          dest="out_dir",
                          type="string",
                          help="Root output directory of cufflinks runs",
                          default=None)
    opt_parser.add_option("--check",
                          dest="check",
                          action="store_true",
                          help="""Will check samples that are not done and print
                                  out which need to still be run""",
                         default=False)
    opt_parser.add_option("--num_processes",
                          dest="num_processes",
                          type="int",
                          help="""If running locally, indicate the number of
                                  processes to batch. Def=%d""" % DEF_NUM_PROCESSES,
                          default=None)
    opt_parser.add_option("--nice",
                          dest="nice",
                          action="store_true",
                          help="If running locally, run using nice",
                          default=False)
    opt_parser.add_option("--LSF",
                          dest="run_lsf",
                          action="store_true",
                          help="Run through LSF",
                          default=None)
    opt_parser.add_option("--print_cmd",
                          dest="print_cmd",
                          action="store_true",
                          help="Print the commands to run, but do not run",
                          default=False)


    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("--out_dir")
    opt_parser.check_required("--txt_ref")

    out_dir = formatDir(options.out_dir)
    if not os.path.exists(out_dir):
        print "Output directory does not exist: %s" % out_dir
        opt_parser.print_help()
        sys.exit(1)

    num_processes = options.num_processes
    run_lsf = options.run_lsf
    nice = options.nice

    quantitate = options.quantitate
    
    print_cmd = options.print_cmd

    force = options.force
    check = options.check
    
    bsub_options = "#!/bin/tcsh\n#BSUB -q week\n#BSUB -R \"rusage[mem=8]\"\n"
    bsub_options += "#BSUB -P cgafolk\n"

    input = open(options.input)

    ctr = 0
    for line in input:
        line = formatLine(line)
        
        s_id, bam = line.split("\t")

        # Make subdir
        subdir = out_dir + "/" + s_id + "_cufflinks"
        if not os.path.exists(subdir):
            os.mkdir(subdir)

        # Check for existence
        file_is_present = False
        try:
            if os.path.getsize(subdir + "/transcripts.gtf") == 0:
                if check:
                    print "Need to run %s" % s_id
            else:
                file_is_present = True
        except: # File doesn't exist
            if check:
                print "Need to run %s" % s_id

        if check:
            continue

        if not force:
            if file_is_present:
                continue
    
        ctr += 1

        cmd = "%s -o %s " % (CUFF_EXEC, subdir)
    
        if quantitate:
            cmd += "-G %s " % options.txt_ref
        else:
            cmd += "-g %s " % options.txt_ref
        cmd += "-u %s" % bam

        if num_processes:
            if nice:
                cmd = "nice " + cmd

            if print_cmd:
                print cmd
                continue

            if ctr % num_processes == 0:
                print cmd
                os.system(cmd)
            else:
                print cmd
                Popen(cmd, shell=True, executable=SHELL)
        else:

            if print_cmd:
                print cmd
                continue

            tmp_file = "%s/tmp_cuff_%s.txt" % (os.curdir,
                                               s_id)
            runLSF(cmd,
                   "%s.cufflinks.bsub.out" % s_id,
                   "cuff_%s" % s_id,
                   "week",
                   tmp_file_name=tmp_file)


    input.close()
			
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
