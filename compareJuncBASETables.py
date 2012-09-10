#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# compareJuncBASETables.py
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

from makeNonRedundantAS import findLargestRegion, updateRedundantDictionary
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
    opt_parser.add_option("--table1",
                          dest="table1",
                          type="string",
                          help="First JuncBASE table",
                          default=None)
    opt_parser.add_option("--table2",
                          dest="table2",
                          type="string",
                          help="""Second JuncBASE table to compare with the first
                                  one""",
                          default=None)
    opt_parser.add_option("--out_prefix",
                          dest="out_prefix",
                          type="string",
                          help="""Prefix to three output tables giving: (1)
                                  events in both tables, (2) events in 1 but not
                                  2, (3) events in 2 but not 1.""",
                          default=None)
#   opt_parser.add_option("--out_1_not_2",
#                         dest="out_1_not_2",
#                         type="string",
#                         help="""JuncBASE table of events in table 1 but not
#                                 table 2.""",
#                         default=None)
#   opt_parser.add_option("--out_2_not_1",
#                         dest="out_2_not_1",
#                         type="string",
#                         help="""JuncBASE table of events in table 2 but not
#                                 table 1.""",
#                         default=None)
#   opt_parser.add_option("--out_1_and_2",
#                         dest="out_1_and_2"
#                         type="string",
#                         help="""JuncBASE table of events in table 2 but not
#                                 table 1.""",
#                         default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--table1")
    opt_parser.check_required("--table2")
    opt_parser.check_required("--out_prefix")


    table1 = open(options.table1)
    table2 = open(options.table2)

    out_1_and_2 = open("%s_1_and_2.txt" % options.out_prefix, "w")
    out_1_not_2 = open("%s_1_not_2.txt" % options.out_prefix, "w")
    out_2_not_1 = open("%s_2_not_1.txt" % options.out_prefix, "w")

    as_type2redundantGroup2event = {}

    # Insert events from table 1
    for event in table1:
        event = formatLine(event)
        
        if event.startswith("#"):
            # Print headers
            out_1_and_2.write(event + "\n")
            out_1_not_2.write(event + "\n")
            out_2_not_1.write(event + "\n")
            continue

        buildDictionary(event, "1")
    table1.close()

    # Now insert events from table 2
    for event in table2:
        event = formatLine(event)

        buildDictionary(event, "2")
    table2.close()

    # Now go through dictionary and output events
    for as_type in as_type2redundantGroup2event:
        for rGroup in as_type2redundantGroup2event[as_type]:
            which_tables = getTableNums(as_type2redundantGroup2event[as_type][rGroup])

            if which_tables == [1,2]:
                this_event = as_type2redundantGroup2event[as_type][rGroup].pop()
                this_event_list = this_event.split("\t")
                out_1_and_2.write("\t".join(this_event_list[:-1]) + "\n")
            elif which_tables = [1]:
                for event in as_type2redundantGroup2event[as_type][rGroup]:
                    this_event_list = this_event.split("\t")
                    out_1_not_2.write("\t".join(this_event_list[:-1]) + "\n")
            elif which_tables = [2]:
                for event in as_type2redundantGroup2event[as_type][rGroup]:
                    this_event_list = this_event.split("\t")
                    out_2_not_1.write("\t".join(this_event_list[:-1]) + "\n")
            else:
                print "Error in redundant group: %s" % rGroup

    out_1_and_2.close()
    out_1_not_2.close()
    out_2_not_1.close()
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildDictionary(event, which_table):
    
    line_list = event.split("\t")

    as_type = line_list[1]

    chr = line_list[3]
    strand = line_list[4]

    # Add on table num to the end fo the event
    event = event + "\t%s" % which_table
    
    if as_type == "mutually_exclusive":
        # Largest region will be in exclusion junctions
        group_start, group_end = findLargestRegion(",".join([line_list[5],
                                                             line_list[6]]))
    else:
        group_start, group_end = findLargestRegion(line_list[5])

    redundantRegion = (chr, strand, group_start, group_end)

    updateRedundantDictionary(as_type2redundantGroup2event, as_type,
                              redundantRegion, event)

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getTableNums(event_set):
    tables = set([])

    for event in event_set:
        event_list = event.split("\t")

        tables.add(int(event_list[-1]))

    table_list = list(tables)

    return table_list.sort()

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
