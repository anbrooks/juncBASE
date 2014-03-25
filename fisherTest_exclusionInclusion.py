#!/lab/64/bin/python
# fisherTest_exclusionInclusion.py
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import pdb

import rpy2.robjects as robjects

r = robjects.r
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
#############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("-f",
                          dest="infile",
                          type="string",
                          help="""Input file where last 4 columns indicate
                                  exclusion and inclusion counts between two
                                  samples.""",
                          default=None)
    opt_parser.add_option("--method",
                          dest="method",
                          type="string",
                          help="""Correction Method: "BH" - Benjamini & Hochberg,
                                  "bonferroni".  Must select these strings as
                                  the option""",
                          default=None)
    opt_parser.add_option("--filter",
                          dest="filter",
                          action="store_true",
                          help="""Filters lines that contain zeros in either
                                  sample 1 or in sample 2.""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-f")
    opt_parser.check_required("--method")
    opt_parser.check_required("--filter")

    file = open(options.infile)

    method = options.method
    if method != "BH" and method != "bonferroni":
        print "Wrong correction method."
        opt_parser.print_help()
        sys.exit(1)

    filter = options.filter

    # List of lines in input file
    infile_lines_preFormat = file.readlines()      

    infile_lines_w_zeros = map(formatLine, infile_lines_preFormat)

    # Remove headerline
    if infile_lines_w_zeros[0].startswith("#"):
        infile_lines_w_zeros.pop(0)

    if filter: 
        infile_lines = removeZeros(infile_lines_w_zeros)
    else:
        infile_lines = infile_lines_w_zeros


    # Parallel list of p_values
    p_vals = [None for i in range(len(infile_lines))]

    # Calculate fisher test for every line and store p_value
    for i in range(len(infile_lines)):
        this_line = infile_lines[i]

        line_list = this_line.split("\t")

        A1 = int(line_list[-4])
        A2 = int(line_list[-3])
        B1 = int(line_list[-2])
        B2 = int(line_list[-1])

        # I sort of get why this indexing works, but the robjects is still
        # unclear to me.  The first [0] indexes to the p_value of the fisher
        # test.  The next[0] then gets the value in the Rvector
        p_vals[i] = robjects.r['fisher.test'](r.matrix(robjects.IntVector([A1,A2, B1, B2]), 
                                              nrow=2))[0][0]

        
    # Get adjusted p_values with Benjamini Hochberg method
    adj_p_vals_rVec = robjects.r['p.adjust'](robjects.FloatVector(p_vals),
                                             method) 

    adj_p_vals = []
    
    for p_val in adj_p_vals_rVec:
        adj_p_vals.append(p_val)

    # Print new results, sorted from lowest p to highest p
    adj_p_vals_sort = list(adj_p_vals)
    adj_p_vals_sort.sort()


    used_line = [False for i in range(len(infile_lines))]

    for this_pval in adj_p_vals_sort:
        for i in range(len(infile_lines)):
            if adj_p_vals[i] == this_pval and used_line[i] == False:
                print infile_lines[i] + "\t" + repr(p_vals[i]) + "\t" +\
                      repr(adj_p_vals[i]) 

                used_line[i] = True

#                del adj_p_vals[i]
#                del p_vals[i]
#                del infile_lines[i]
    
                break
			
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

def removeZeros(lines_w_zeros):
    return_lines = []

    for line in lines_w_zeros:
        line_list = line.split("\t")

        A1 = int(line_list[-4])
        A2 = int(line_list[-3])
        B1 = int(line_list[-2])
        B2 = int(line_list[-1])

        # If one sample has zeroes in both columns, then do not consider this
        # line.

        if A1 == 0 and A2 == 0:
            continue
        if B1 == 0 and B2 == 0:
            continue
        
        return_lines.append(line)

    return return_lines

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
