#!/lab/64/bin/python
# compareSampleSets_getASEvents.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Will do test to compare two sets of groups in the data.

   Input will be an output file from clusterASExons2.py
"""

import sys
import optparse 
import pdb

import rpy2.robjects as robjects
from helperFunctions import updateDictOfLists

r = robjects.r
# Suppresses warnings
robjects.r["options"](warn=-1)
#############
# CONSTANTS #
#############
NA = "NA"
SIGN_CUTOFF = 0.05
DEF_THRESH = 25
DEF_DPSI_THRESH = 5.0

#SAMP_SET_THRESH = 3

DEF_TEST = "Wilcoxon"
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
                          help="Resulting file from clusterASExons2.py",
                          default=None)
    opt_parser.add_option("--left_intron",
                          dest="left_input",
                          type="string",
                          help="""Resulting file from clusterASExons2.py, which
                                  contains the exclusion and inclusion counts
                                  for just the left side of an intron retention
                                  event.""",
                          default=None)
    opt_parser.add_option("--right_intron",
                          dest="right_input",
                          type="string",
                          help="""Resulting file from clusterASExons2.py, which
                                  contains the exclusion and inclusion counts
                                  for just the right side of an intron retention
                                  event.""",
                          default=None)
    opt_parser.add_option("--all_psi_output",
                          dest="all_psi_output",
                          type="string",
                          help="""Output file that will contain the PSI values
                                  for all events and samples. The last two
                                  columns will correspond to the raw-pvalue and
                                  corrected p-value.""",
                          default=None)
    opt_parser.add_option("--thresh",
                          dest="threshold",
                          type="int",
                          help="""Threshold for minimum number of total reads
                                  in an event. Default=%d""" % DEF_THRESH,
                          default=DEF_THRESH)
    opt_parser.add_option("--mt_correction",
                          dest="mt_method",
                          type="string",
                          help="""Multiple testing correction Method: "BH" - Benjamini & Hochberg,
                                  "bonferroni".  Must select these strings as
                                  the option""",
                          default=None)
    opt_parser.add_option("--which_test",
                          dest="which_test",
                          type="string",
                          help="""Which test to use. Either "t-test" or
                                  "Wilcoxon". Default=%s""" % DEF_TEST,
                          default=DEF_TEST)
    opt_parser.add_option("--as_dPSI_thresh",
                          dest="as_dPSI_thresh",
                          type="float",
                          help="""Minimum PSI difference between the maximum
                                  and minimum PSI values for a given event to be
                                  called an alternatively spliced event. This
                                  should probably be less than the dPSI
                                  threshold used to filter significantly
                                  associated events. Default=%s""" % DEF_DPSI_THRESH,
                          default=DEF_DPSI_THRESH)
    opt_parser.add_option("--sample_set1",
                          dest="sample_set1",
                          type="string",
                          help="""Comma delimited list of samples in set 1.
                                  Names must be in header columns of input
                                  files.""",
                          default=None)
    opt_parser.add_option("--sample_set2",
                          dest="sample_set2",
                          type="string",
                          help="""Comma delimited list of samples in set 2.
                                  Names must be in header columns of input
                                  files.""",
                          default=None)
#   opt_parser.add_option("--html_out_dir",
#                         dest="html_out_dir",
#                         type="string",
#                         help="Directory to output html report.",
#                         default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("--all_psi_output")
    opt_parser.check_required("--mt_correction")
    opt_parser.check_required("--sample_set1")
    opt_parser.check_required("--sample_set2")

    input_file = open(options.input_file)
    left_input_file_name = options.left_input
    right_input_file_name = options.right_input
    sum_thresh = options.threshold

    as_dPSI_thresh = options.as_dPSI_thresh

    left_input_file = None
    right_input_file = None
    if left_input_file_name is None:
        print "Warning: No intron retention file given as input.  Will not calculate IR events."
    else:
        left_input_file = open(left_input_file_name)
        right_input_file = open(right_input_file_name)
    
    all_psi_output = open(options.all_psi_output, "w")

    html_out_dir = options.html_out_dir
    if html_out_dir:
        html_out_dir = formatDir(html_out_dir)

    method = options.mt_method
    if method != "BH" and method != "bonferroni":
        print "Wrong method indicated."
        opt_parser.print_help()
        sys.exit(1)

    which_test = options.which_test 
    if which_test != "Wilcoxon" and which_test != "t-test":
        print "Wrong method indicated."
        opt_parser.print_help()
        sys.exit(1)

    if which_test == "Wilcoxon":
        which_test = "wilcox.test"
    if which_test == "t-test":
        which_test = "t.test"

    sample_set1 = options.sample_set1.split(",")
    sample_set2 = options.sample_set2.split(",")

    # The threshold for the number of samples that need to have expressed AS
    # events in order to consider testing
    samp_set_thresh1 = float(len(sample_set1)) / 2
    samp_set_thresh2 = float(len(sample_set2)) / 2

    idx2sample = {}

    # {event_type:[pval]}
    event_type2pvals = {}

    # {event_type:(set1_medianPSI, set2medianPSI),]}
    event_type2PSI_vals_4_set = {}

    # {event::pval_idx}
    event2idx = {}

    # {event:{col:psi}}
    event2col2psi = {}

    # {event:{col:sum_counts}}
    event2col2sum = {}

    header = None
    total_samples = None
    for line in input_file:
        line = formatLine(line)

        if line.startswith("#"):
            header = line
            headerList = header.split("\t")
            sampleList = headerList[11:]
            # Get sample idx
            for i in range(len(sampleList)):
                idx2sample[i] = sampleList[i]
#           for sample in sample_set1:
#               idx2sample[sampleList.index(sample)] = sample
#           for sample in sample_set2:
#               idx2sample[sampleList.index(sample)] = sample
            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:11])
        counts = line_list[11:]

        event_type = getEventType(event)
        if event_type not in event_type2pvals:
            event_type2pvals[event_type] = []
            event_type2PSI_vals_4_set[event_type] = []
        total_samples = len(counts)

        # Fill PSI dict
        min_psi = 200
        max_psi = -1
        for i in range(len(counts)):
            (psi, sum_ct) = getPSI_sample_sum(counts[i], sum_thresh)
            if psi != NA:
                if psi < min_psi:
                    min_psi = psi
                if psi < max_psi:
                    max_psi = psi
            if event in event2col2psi:
                event2col2psi[event][i] = psi
                event2col2sum[event][i] = sum_ct
            else:
                event2col2psi[event] = {i:psi}
                event2col2sum[event] = {i:sum_ct}


        # Compare samples groups together in a wilcoxon rank sum test
        set1_psis = []        
        set2_psis = []
        for j in range(len(counts)):
            [col_excl, col_incl] = map(int,counts[j].split(";"))

            # Both samples have to be non-zero
            if belowThreshold(sum_thresh, col_excl, col_incl):
                continue

            if idx2sample[j] in sample_set1:
                if event2col2psi[event][j] != NA:
                    set1_psis.append(event2col2psi[event][j])
            elif idx2sample[j] in sample_set2:
                if event2col2psi[event][j] != NA:
                    set2_psis.append(event2col2psi[event][j])

        if len(set1_psis) <= samp_set_thresh1 or len(set2_psis) <= samp_set_thresh2:
            continue

        if (max_psi - min_psi) < as_dPSI_thresh:
            continue
        
        event_type2PSI_vals_4_set[event_type].append((robjects.r['median'](robjects.FloatVector(set1_psis)),
                                                      robjects.r['median'](robjects.FloatVector(set2_psis))))

        # Calculate p-val for intron retention later
        if event_type == "intron_retention":
            continue

        cur_len = len(event_type2pvals[event_type])
#        cur_len2 = len(event_type2col2pvals[event_type][j])

#           if event in event2pairs2idx:
#               event2pairs2idx[event][(0,j)] = cur_len
#           else:
#               event2pairs2idx[event] = {(0,j):cur_len}	

#           if event in event2col2idx:
#               event2col2idx[event][j] = cur_len2
#           else:
#               event2col2idx[event] = {j:cur_len2}
#         
       
        try: 
            raw_pval = robjects.r[which_test](robjects.FloatVector(set1_psis),
                                          robjects.FloatVector(set2_psis))[2][0]
        except:
            continue

        if robjects.r["is.nan"](raw_pval)[0]:
            continue

        event_type2pvals[event_type].append(raw_pval)
        event2idx[event] = cur_len

    # Now calculate intron retention
    if left_input_file:
        left_events2counts = getIntronLeftRightCounts(left_input_file)
        right_events2counts = getIntronLeftRightCounts(right_input_file)
    else:
        left_events2counts = {}
        right_events2counts = {}

    for event in left_events2counts:
        if event not in right_events2counts:
            continue

        set1_psis_left = []        
        set2_psis_left = []
        set1_psis_right = []        
        set2_psis_right = []

        left_min_psi = 200
        left_max_psi = -1
        right_min_psi = 200
        right_max_psi = -1
        for j in range(total_samples):
            [left_col_excl, left_col_incl] = map(int,left_events2counts[event][j].split(";"))
            [right_col_excl, right_col_incl] = map(int,right_events2counts[event][j].split(";"))

            # Both samples have to be non-zero
            if (belowThreshold(sum_thresh, left_col_excl, left_col_incl)
                               or
                belowThreshold(sum_thresh, right_col_excl, right_col_incl)):
                continue

            (left_psi, sum_ct) = getPSI_sample_sum(left_events2counts[event][j], sum_thresh)
            (right_psi, sum_ct) = getPSI_sample_sum(right_events2counts[event][j], sum_thresh)

            if left_psi != NA:
                if left_psi < left_min_psi:
                    left_min_psi = left_psi
                if left_psi > left_max_psi:
                    left_max_psi = left_psi

            if right_psi != NA:
                if right_psi < right_min_psi:
                    right_min_psi = right_psi
                if right_psi > right_max_psi:
                    right_max_psi = right_psi

            if idx2sample[j] in sample_set1:
                if left_psi != NA:
                    set1_psis_left.append(left_psi)
                if right_psi != NA:
                    set1_psis_right.append(right_psi)
            elif idx2sample[j] in sample_set2:
                if left_psi != NA:
                    set2_psis_left.append(left_psi)
                if right_psi != NA:
                    set2_psis_right.append(right_psi)

        if len(set1_psis_left) <= samp_set_thresh1 or len(set1_psis_right) <= samp_set_thresh1\
            or len(set2_psis_left) <= samp_set_thresh2 or len(set2_psis_right) <= samp_set_thresh2:
            continue

        if (left_max_psi - left_min_psi) < as_dPSI_thresh:
            continue
        if (right_max_psi - right_min_psi) < as_dPSI_thresh:
            continue

        cur_len = len(event_type2pvals["intron_retention"])

        try:
            left_pval = robjects.r[which_test](robjects.FloatVector(set1_psis_left),
                                               robjects.FloatVector(set2_psis_left))[2][0]
                    
            right_pval = robjects.r[which_test](robjects.FloatVector(set1_psis_right),
                                                robjects.FloatVector(set2_psis_right))[2][0]
        except:
            continue

        if robjects.r["is.nan"](left_pval)[0] or robjects.r["is.nan"](right_pval)[0]:
            continue
        else:
            combined_pval = (left_pval + right_pval) - left_pval * right_pval

        event_type2pvals["intron_retention"].append(combined_pval)
        event2idx[event] = cur_len

    # All pairs have been evaluated, so now do multiple testing correction on
    # everything
    event_type2adjusted_pvals = {}
    event_type2col2adjusted_pvals = {}

    for event_type in event_type2pvals:
        event_type2adjusted_pvals[event_type] = robjects.r['p.adjust'](robjects.FloatVector(event_type2pvals[event_type]),
                                                                       method) 
    
    # Now go through all events and print out pvals
    all_psi_output.write(header + "set1_med_psi\tset2_med_psi\tdeltaPSI\traw_pval\tcorrected_pval\n")

    for event in event2idx:
        event_type = getEventType(event)

        this_idx = event2idx[event]
        if this_idx == NA:
            psi_vals = []
            for i in range(total_samples):
                psi_vals.append(event2col2psi[event][i])

            outline = "%s\t%s\tNA\tNA\n" % (event, 
                                            "\t".join(psi_vals))

            all_psi_output.write(outline)
            continue

        psi_vals = []
        for i in range(total_samples):
            psi_vals.append(event2col2psi[event][i])

        outline = "%s\t%s" % (event, 
                                "\t".join(psi_vals))

        # Add median PSI and delta PSI values
        outline += "\t%.2f\t%.2f\t%.2f" % (event_type2PSI_vals_4_set[event_type][this_idx][0],
                                           event_type2PSI_vals_4_set[event_type][this_idx][1],
                                           event_type2PSI_vals_4_set[event_type][this_idx][1] -
                                           event_type2PSI_vals_4_set[event_type][this_idx][0])

        outline += "\t%f\t%f\n" % (event_type2pvals[event_type][this_idx],
                                   event_type2adjusted_pvals[event_type][this_idx])

        all_psi_output.write(outline)

    all_psi_output.close()

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############

def belowThreshold(sum_thresh, col_excl, col_incl):

    if col_excl + col_incl < sum_thresh:
        return True

    return False

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir
    
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getPSI_sample_sum(excl_incl_ct_str, sum_thresh):
    if excl_incl_ct_str == NA:
        return NA, NA

    excl_str, incl_str = excl_incl_ct_str.split(";")

    try:
        excl = float(excl_str)
        incl = float(incl_str)
    except:
        print "Warning: Bad counts for PSI"
        return (NA, NA)

    if (excl + incl) <  sum_thresh:  
        return (NA, NA)
   
    psi = (incl/(incl + excl)) * 100

    psi_str = "%.2f" % psi
    total_ct = "%d" % int(excl + incl)

    return (psi_str, total_ct)

     
def getEventType(event):
    return event.split("\t")[1]
    
def getIntronLeftRightCounts(file):

    intron_event2counts = {}
    
    for line in file:
        line = formatLine(line)

        if line.startswith("#"):
            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:11])
        counts = line_list[11:]

        # If the reference is NA, then do not calculate
        if counts[0] == NA:
            continue

        intron_event2counts[event] = counts

    return intron_event2counts
        
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
