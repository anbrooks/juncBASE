#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# getSimpleAS.py
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
SIMPLE_EVENTS = set(["cassette",
                 "alternative_donor",
                 "alternative_acceptor"])
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
                          dest="juncBASE_table",
                          type="string",
                          help="JuncBASE table",
                          default=None)
    opt_parser.add_option("-o",
                          dest="out_table",
                          type="string",
                          help="""juncBASE_table only containing simple cassette exon
                                  events (no additional alternative splice site
                                  usage) and alternative 5' and alternative 3'
                                  splice sites with two choices""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")

    in_table = open(options.juncBASE_table)
    out_table = open(options.out_table, "w")

    for line in in_table:
        line = formatLine(line)

        # Header line
        if line.startswith("#"):
            out_table.write(line + "\n") 
            continue

        lineList = line.split("\t")
    
        if lineList[1] not in SIMPLE_EVENTS:
            continue

        if lineList[1] == "cassette":
            # Only keep events with one skipping junction
            if ";" in lineList[5]:
                continue

        if lineList[1] == "alternative_donor" or lineList[1] == "alternative_acceptor":
            if ";" in lineList[5]:
                continue
            if "," in lineList[5]:
                continue
        
        out_table.write(line + "\n") 
			
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
