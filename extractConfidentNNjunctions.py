#!/lab/64/bin/python
# extractConfidentNNjunctions.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Will get only the NN titles that contain the NN junctions that were
considered confident.
"""

import sys
import optparse 
import os

from subprocess import Popen
from helperFunctions import runCmd
#############
# CONSTANTS #
#############
SCRIPT = os.path.realpath("./getSpecificSequences.py")
if not os.path.exists(SCRIPT):
    print "ERROR: getSpecificSequences.py needs to be in the same directory."
    sys.exit(1)
SHELL = "/bin/tcsh"

NUM_PROC = 20
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
                          dest="split_nn_dir",
                          type="string",
                          help="""Directory containing the NN junctions split
                                  into multiple fasta files.""",
                          default=None)
    opt_parser.add_option("-t",
                          dest="titles",
                          type="string",
                          help="""File containing title names of all junctions
                                  that will be extracted.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-d")
    opt_parser.check_required("-t")

    file_counter = 0

    nn_dir = options.split_nn_dir

    if not nn_dir.endswith("/"):
        nn_dir += "/"

    title_file = options.titles

    # Remove all tmp files currently in the directory
    cmd = "rm tmp*.txt"
    runCmd(cmd, SHELL)

    for fa_file in os.listdir(nn_dir):
        if not fa_file.endswith(".fa"):
            continue

        file_counter += 1
        tmp_file = "tmp%d.txt" % file_counter
           
        cmd = "python %s -f %s -t %s > %s" % (SCRIPT,
                                              nn_dir + fa_file,
                                              title_file,
                                              tmp_file)

        if file_counter % NUM_PROC == 0:
            runCmd(cmd, SHELL)
        else:
            p = Popen(cmd, shell=True, executable=SHELL)
   
     
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
