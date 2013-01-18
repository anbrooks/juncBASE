#!/lab/64/bin/python
# getSpecificSequences.py
# Author: Angela Brooks
# Program Completion Date:
# Description: Take a list of titles and a fasta file and only outputs
# sequences that are in the fasta file.
# Modification Date(s):

import sys
import optparse 
import pdb

from Bio import SeqIO
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
    opt_parser.add_option("-f",
                          dest="fasta_file",
                          type="string",
                          help="Fasta file of sequences.",
                          default=None)
    opt_parser.add_option("-t",
                          dest="titles",
                          type="string",
                          help="Titles of sequences that you want to extract.",
                          default=None)
    opt_parser.add_option("-c",
                          dest="isContained",
                          action="store_true",
                          help="""Title input just has to be contained in fasta
                                  record.  Does not have to be an exact
                                  match.""",
                          default=False)
    

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-f")
    opt_parser.check_required("-t")
	
    fasta_file = open(options.fasta_file)
    title_file = open(options.titles)

    isContained = options.isContained
	
    # Put titles into a dictionary for searching
    # Title maps to arbitrary number, 0
    title_set = set([])

    for line in title_file:
        line = formatLine(line)

        if line.startswith(">"):
            line = line[1:] # Chomp off the >

        title_set.add(line)

    # Check for title in dictionary
    for record in SeqIO.parse(fasta_file, "fasta"):

        title = record.id

        if isContained:
            for t in title_set:
                if t in title:
                    print record.format("fasta"),
                    break
        else:
            if title in title_set:
                print record.format("fasta"),
	
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
