#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# run_preProcess_step3_by_chr.py
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
from helperFunctions import runCmd, runLSF
#############
# CONSTANTS #
#############
SHELL = "/bin/tcsh"

BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/preProcess_getASEventReadCounts_step3.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: preProcess_getASEventReadCounts_step3.py needs to be in the same directory."
    sys.exit(1)

DEF_OVERHANG = 6

DEF_NUM_PROCESSES = 2
DEF_GROUP = "cgafolk"
DEF_QUEUE = "week"
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
    opt_parser.add_option("--input_dir",
                          dest="input_dir",
                          type="string",
                          help="Root directory of all input files.",
                          default=None)
    opt_parser.add_option("--LSF",
                          dest="lsf_group",
                          type="string",
                          help="""Launches jobs on LSF by default. It will
                                  use this specified group. DEF=%s""" % DEF_GROUP,
                          default=DEF_GROUP)
    opt_parser.add_option("--lsf_queue",
                          dest="lsf_queue",
                          type="string",
                          help="""If launching jobs on LSF, it will use the
                                  specified queue. DEF=%s""" % DEF_QUEUE,
                          default=DEF_QUEUE)
    opt_parser.add_option("--num_processes",
                          dest="num_processes",
                          type="int",
                          help="""If jobs should be run locally, indicate 
                                  by putting the number of processes to batch.""",
                          default=None)
    opt_parser.add_option("--min_overhang",
                          dest="min_overhang",
                          type="int",
                          help="""Minimum overhang used to determine a junction
                                  alignment (Used as input to find intron-exon
                                  junctions. Default=%d""" % DEF_OVERHANG,
                         default=DEF_OVERHANG)
    opt_parser.add_option("--force",
                          dest="force",
                          action="store_true",
                          help="""By default, will check for the existence of the output
                                  file. If it is non-zero, then it will not run
                                  the sample. This option forces runs of
                                  everything""",
                         default=False)
    opt_parser.add_option("--check",
                          dest="check",
                          action="store_true",
                          help="""Will check samples that are not done and print
                                  out which need to still be run""",
                         default=False)
    opt_parser.add_option("--nice",
                          dest="nice",
                          action="store_true",
                          help="""Will run locally, using nice""",
                         default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--input_dir")

    input_dir = formatDir(options.input_dir)

    lsf_group = options.lsf_group
    lsf_queue = options.lsf_queue
    num_processes = options.num_processes

    force = options.force
    check = options.check
    nice = options.nice

    # Will use the tmp files in the input directory to determine chromosomes to
    # process.
    chr_list = []
    for bed_file in os.listdir(input_dir):
        if os.path.isdir(input_dir + "/" + bed_file):
            continue

        if not bed_file.endswith(".bed"):
            continue

        first_split = bed_file.split("preProcess")[0]
        second_split = first_split.split("tmp")[-1]
        chr_list.append(second_split.strip("_"))
    
    # Now go through each sample directory to get to subdirectory, then run
    # command
    ctr = 0 # For num processes
    for sample_dir in os.listdir(input_dir):
        if not os.path.isdir(input_dir + "/" + sample_dir):
            continue

        for this_chr in chr_list:
            expected_out_file = "%s/%s/%s_%s/%s_%s_intron_exon_junction_counts.txt" % (input_dir,
                                                                                       sample_dir,
                                                                                       sample_dir, this_chr,
                                                                                       sample_dir, this_chr)

            file_is_present = False
            try:
                if os.path.getsize(expected_out_file) == 0:
                    if check:
                        print "File is empty: %s" % expected_out_file 
                else:
                    file_is_present = True
            except:
                if check:
                    print "Does not exist: %s" % expected_out_file
           
            if check: 
                continue

            if not force:
                if file_is_present:
                    continue 

            ctr += 1

            cmd = "python %s " % SCRIPT
            cmd += "-i %s/%s " % (input_dir, sample_dir)
            cmd += "-n %s_%s " % (sample_dir, this_chr)
            cmd += "-t %s/tmp_%s_preProcess_getASEventReadCounts_step2.bed " % (input_dir,
                                                                                this_chr)
            cmd += "--min_overhang %d" % options.min_overhang

            if num_processes:
                if nice:
                    cmd = "nice " + cmd
                if ctr % num_processes == 0:
                    os.system(cmd)
                else:
                    print cmd
                    Popen(cmd, shell=True, executable=SHELL)
            else:
                tmp_file = "%s/tmp.txt" % os.curdir
                runLSF(cmd,  
                       "%s_%s.preProcess_3.bsub.out" % (sample_dir, this_chr),
                       "%s_%s" % (sample_dir, this_chr),
                       lsf_queue,
                       group=lsf_group,
                       tmp_file_name=tmp_file)

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
