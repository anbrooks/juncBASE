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
    opt_parser.add_option("--compare_pval",
                          dest="compare_pval",
                          action="store_true",
                          help="""When the same event has been identified, the
                                  p-value of the significance for table 1 and
                                  table 2 will be given in the last two
                                  columns.""",
                          default=False)
    opt_parser.add_option("--subset",
                          dest="subset_prefix",
                          type="string",
                          help="""Optional: Prefix of a table that will contain a subset
                                  of table1 with only events found in table2""",
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

    subset = options.subset_prefix
    compare_pval = options.compare_pval    

    subset_file = None
    if subset:
        subset_file = open(subset + "_subset.txt", "w")

    table1 = open(options.table1)
    table2 = open(options.table2)

    out_1_and_2 = open("%s_1_and_2.txt" % options.out_prefix, "w")
    out_1_not_2 = open("%s_1_not_2.txt" % options.out_prefix, "w")
    out_2_not_1 = open("%s_2_not_1.txt" % options.out_prefix, "w")

    as_type2chr2strand2redundantGroup2event = {}

    # Insert events from table 1
    for event in table1:
        event = formatLine(event)
        
        if event.startswith("#"):
            # Print headers
            out_1_and_2.write(event + "\n")
            out_1_not_2.write(event + "\n")
            out_2_not_1.write(event + "\n")
            continue

        buildDictionary(as_type2chr2strand2redundantGroup2event, event, "1")
    table1.close()

    table2_events = set([])
    # Now insert events from table 2
    for event in table2:
        if event.startswith("#"):
            continue    

        event = formatLine(event)

        buildDictionary(as_type2chr2strand2redundantGroup2event, event, "2")

        if subset_file:
            lineList = event.split("\t")
            table2_events.add((lineList[1],
                               lineList[3],
                               lineList[4],
                               lineList[5],
                               lineList[6]))

    table2.close()

    # Now go through dictionary and output events
    for as_type in as_type2chr2strand2redundantGroup2event:
        for chr in as_type2chr2strand2redundantGroup2event[as_type]:
            for strand in as_type2chr2strand2redundantGroup2event[as_type][chr]:
                for rGroup in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                    which_tables, pvals = getTableNums(as_type2chr2strand2redundantGroup2event[as_type][chr][strand][rGroup],
                                                       compare_pval)
                    if which_tables == [1,2]:
                        this_event = as_type2chr2strand2redundantGroup2event[as_type][chr][strand][rGroup].pop()
                        this_event_list = this_event.split("\t")
                        # Only the event information will be given as output 
                        out_1_and_2.write("\t".join(this_event_list[:11])) 
                        if compare_pval:
                            out_1_and_2.write("\t%.3f\t%.3f" % (pvals[0],
                                                                pvals[1]))
                        out_1_and_2.write("\n")
                    elif which_tables == [1]:
                        for this_event in as_type2chr2strand2redundantGroup2event[as_type][chr][strand][rGroup]:
                            this_event_list = this_event.split("\t")
                            out_1_not_2.write("\t".join(this_event_list[:-1]) + "\n")
                    elif which_tables == [2]:
                        for this_event in as_type2chr2strand2redundantGroup2event[as_type][chr][strand][rGroup]:
                            this_event_list = this_event.split("\t")
                            out_2_not_1.write("\t".join(this_event_list[:-1]) + "\n")
                    else:
                        print "Error in redundant group: %s" % repr(rGroup)
                

    out_1_and_2.close()
    out_1_not_2.close()
    out_2_not_1.close()

    # Will print subset if option is given
    if subset_file:
        table1 = open(options.table1)
        for line in table1:
            if line.startswith("#"):
                subset_file.write(line)
                continue

            line = formatLine(line)
            lineList = line.split("\t")
            if (lineList[1],
                lineList[3],
                lineList[4],
                lineList[5],
                lineList[6]) in table2_events:
                subset_file.write(line + "\n")

        table1.close()

    if not subset is None:
        subset_file.close()
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildDictionary(as_type2chr2strand2redundantGroup2event, event, which_table):
   
    redundantRegion, as_type = get_rGroup_as_event(event)

    # Add on table num to the end fo the event
    event = event + "\t%s" % which_table

    updateRedundantDictionary(as_type2chr2strand2redundantGroup2event, as_type,
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

def get_rGroup_as_event(event):
    line_list = event.split("\t")

    as_type = line_list[1]

    chr = line_list[3]
    strand = line_list[4]

    if as_type == "mutually_exclusive":
        # Largest region will be in exclusion junctions
        group_start, group_end = findLargestRegion(",".join([line_list[5],
                                                             line_list[6]]))
    else:
        group_start, group_end = findLargestRegion(line_list[5])

    redundantRegion = (chr, strand, group_start, group_end)

    return redundantRegion, as_type


def getTableNums(event_set, compare_pval=False):
    tables = set([])
    pvals = None

    if compare_pval:
        tab1_pval = 1.0
        tab2_pval = 1.0

    for event in event_set:
        event_list = event.split("\t")
        table_number = int(event_list[-1])

        tables.add(table_number)

        if compare_pval:
            pval = float(event_list[-3])

            if table_number == 1:
                if pval < tab1_pval:
                    tab1_pval = pval
            else:
                if pval < tab2_pval:
                    tab2_pval = pval

    table_list = list(tables)
    table_list.sort()

    if compare_pval:
        pvals = [tab1_pval, tab2_pval]

    return table_list, pvals

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
