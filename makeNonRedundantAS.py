#!/lab/64/bin/python
# makeNonRedundantAS.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
"""From a table containing all AS events with p_vals, removes redundant events.
Assumes that all events were already considered significant.
"""

import sys
import optparse 
import pdb

from helperFunctions import coordsOverlap
from getASEventReadCounts import convertCoordStr

import rpy2.robjects as robjects
#############
# CONSTANTS #
#############
INFINITY = 100000000000000000000

NA = "NA"
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
    opt_parser.add_option("--input",
                          dest="input",
                          type="string",
                          help="""Table containing PSI values for events.  Last
                                  column indicates the combined p-value.""",
                          default=None)
    opt_parser.add_option("--output",
                          dest="output",
                          type="string",
                          help="""Resulting table with all redundancies
                                  removed.""",
                          default=None)
    opt_parser.add_option("--use_mad",
                          dest="use_mad",
                          action="store_true",
                          help="""Instead of a p-value, redundancies will be
                                  resolved by taking the event with the largest 
                                  median absolute deviation (MAD). All value
                                  columns should be PSI values.""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--input")
    opt_parser.check_required("--output")

    input_file = open(options.input)
    output_file = open(options.output, "w")

    use_mad = options.use_mad

    as_type2event2pval = {}

    as_type2chr2strand2redundantGroup2event = {}

    for event in input_file:
        event = formatLine(event)
    
        if event.startswith("#"):
            # Print header
            output_file.write(event + "\n")
            continue
        
        buildDictionaries(event, as_type2event2pval,
                          as_type2chr2strand2redundantGroup2event, use_mad)

    input_file.close()

    # Merge redundant events
    changeOccurred = True
    while changeOccurred:
        changeOccurred = mergeRedundantEvents(as_type2chr2strand2redundantGroup2event)

    # Now remove redundant events from as_type2event2pval dictionary
    removeRedundantEvents(as_type2chr2strand2redundantGroup2event, as_type2event2pval)

    # Now print out the remaining events
    for as_type in as_type2event2pval:
        for event in as_type2event2pval[as_type]:
            output_file.write(event + "\n")

    output_file.close()
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildDictionaries(event, as_type2event2pval, as_type2chr2strand2redundantGroup2event,
                      use_mad):

    line_list = event.split("\t")

    if use_mad:
        pval = get_mad(line_list)
    else:
        pval = float(line_list[-1])
    as_type = line_list[1]

    if as_type in as_type2event2pval:
        if event in as_type2event2pval[as_type]:
            print "WARNING: Duplicate event exists."
            if pval > as_type2event2pval[as_type][event]:
                # Lowest pval is selected
                pval = as_type2event2pval[as_type][event]
        as_type2event2pval[as_type][event] = pval
    else:
        as_type2event2pval[as_type] = {event:pval}

    chr = line_list[3]
    strand = line_list[4]

    if as_type == "mutually_exclusive":
        # Largest region will be in exclusion junctions
        group_start, group_end = findLargestRegion(",".join([line_list[5],
                                                             line_list[6]]))    
    else:
        group_start, group_end = findLargestRegion(line_list[5])    

    redundantRegion = (chr, strand, group_start, group_end)

    updateRedundantDictionary(as_type2chr2strand2redundantGroup2event, as_type, redundantRegion,
                              event)

#ef convertCoordStr(coord_str):
#   chr, start_str, end_str = coord_str.split("_")

#   return chr, int(start_str), int(end_str)

def findLargestRegion(coords_string):

    first_list = coords_string.split(";")
    full_list = []
    for elem1 in first_list:
        for elem2 in elem1.split(","):
            full_list.append(elem2)

    leftmost = INFINITY
    rightmost = -1 

    for coord in full_list:
        chr, start, end = convertCoordStr(coord)

        if start < leftmost:
            leftmost = start

        if end > rightmost:
            rightmost = end

    return leftmost, rightmost 

def findMostSignEvent(event2pval, events):
    most_sign_pval = 2
    most_sign_event = None

    for event in events:
        if event2pval[event] < most_sign_pval:
            most_sign_pval = event2pval[event]
            most_sign_event = event

    return most_sign_event

    
def fixMissingVals(vals_str):

    return_vals = []

    for val in vals_str:
        if val != NA:
            return_vals.append(val)

    return return_vals

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def get_mad(line_list):
    vals_str = line_list[11:]
    
    vals_str = fixMissingVals(vals_str)

    vals = map(float, vals_str) 

    abs_mad_val = abs(robjects.r['mad'](robjects.FloatVector(vals))[0])

    # Since the lowest p-value is supposed to be used for selecting the
    # non-redundant set of events, the negative of the larget MAD value will
    # give the event with the largets MAD value.
    return -abs_mad_val

def mergeRedundantEvents(as_type2chr2strand2redundantGroup2event):
    """
    Because of the construction of the dictionary, some of the redundant
    groups are overlapping.
    """
    for as_type in as_type2chr2strand2redundantGroup2event:
        for chr in as_type2chr2strand2redundantGroup2event[as_type]:
            for strand in as_type2chr2strand2redundantGroup2event[as_type][chr]:
                for first_group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                    for second_group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                        if first_group == second_group:
                            continue
                        if coordsOverlap(first_group[2], first_group[3],
                                         second_group[2], second_group[3]):

                            # Add sets together
                            new_set = as_type2chr2strand2redundantGroup2event[as_type][chr][strand][first_group].union(as_type2chr2strand2redundantGroup2event[as_type][chr][strand][second_group])
                            new_start = min(first_group[2], second_group[2])
                            new_end = max(first_group[3], second_group[3])

                            new_region = (first_group[0], first_group[1],
                                          new_start, new_end)

                            if new_region in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                                as_type2chr2strand2redundantGroup2event[as_type][chr][strand][new_region].update(new_set)
                            else:
                                as_type2chr2strand2redundantGroup2event[as_type][chr][strand][new_region] = new_set

                            # Remove old sets and return with flag
                            if new_region != first_group:
                                del as_type2chr2strand2redundantGroup2event[as_type][chr][strand][first_group]
                            if new_region != second_group:
                                del as_type2chr2strand2redundantGroup2event[as_type][chr][strand][second_group]

                            return True

    # If not returned from within the loop, then no change occurred
    return False

def removeRedundantEvents(as_type2chr2strand2redundantGroup2event, as_type2event2pval):
    """
    Returns the as_type2event2pval dictionary with redundant values removed
    """
    for as_type in as_type2chr2strand2redundantGroup2event:
        for chr in as_type2chr2strand2redundantGroup2event[as_type]:
            for strand in as_type2chr2strand2redundantGroup2event[as_type][chr]:
                for redundant_group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                    if len(as_type2chr2strand2redundantGroup2event[as_type][chr][strand][redundant_group]) > 1:
                        most_sign_event = findMostSignEvent(as_type2event2pval[as_type],
                                                            as_type2chr2strand2redundantGroup2event[as_type][chr][strand][redundant_group])
                
                        for event in as_type2chr2strand2redundantGroup2event[as_type][chr][strand][redundant_group]:
                            if event == most_sign_event:
                                continue
                            # Delete all non-significant events
                            del as_type2event2pval[as_type][event] 

def updateRedundantDictionary(as_type2chr2strand2redundantGroup2event, as_type, redundantRegion,
                              event):

    this_chr = redundantRegion[0]
    this_strand = redundantRegion[1]    

    if as_type in as_type2chr2strand2redundantGroup2event:
        if not this_chr in as_type2chr2strand2redundantGroup2event[as_type]:
            as_type2chr2strand2redundantGroup2event[as_type][this_chr] = {this_strand:{redundantRegion: set([event])}}
            return
        if not this_strand in as_type2chr2strand2redundantGroup2event[as_type][this_chr]:
            as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand] = {redundantRegion: set([event])}
            return

        foundOverlap = False
        for redundantGroup in as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand]:
            if coordsOverlap(redundantGroup[2], redundantGroup[3],
                             redundantRegion[2], redundantRegion[3]):
                foundOverlap = True
                # Pick longest region as the key
                region_start = min(redundantGroup[2], redundantRegion[2])
                region_end = max(redundantGroup[3], redundantRegion[3])

                # Copy the existing set of exon_names
                cur_set = set(as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantGroup]) 
                # Add the current event
                cur_set.add(event)

                del as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantGroup]

                # Check for existing new group
                new_group = (redundantGroup[0], redundantGroup[1],
                             region_start, region_end)
                if new_group in as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand]:
                    as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][new_group].update(cur_set)
                else:
                    as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][new_group] = cur_set
                break
        if not foundOverlap:
            # Add this non overlapping region on its own
            as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantRegion] = set([event]) 
    else:
        as_type2chr2strand2redundantGroup2event[as_type] = {this_chr:{this_strand:{redundantRegion: set([event])}}}

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
