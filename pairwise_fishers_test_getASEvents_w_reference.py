#!/lab/64/bin/python
# pairwise_fishers_test_getASEvents.py 
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Will do fisher's test for all pairwise comparsions.

   Input will be an output file from clusterASExons2.py
"""

import sys
import optparse 
import pdb
import os

import rpy2.robjects as robjects
from helperFunctions import updateDictOfLists
from getASEventReadCounts import convertCoordStr, getAD_AA_isoform_lengths, determineAltStartOrEnd
r = robjects.r
#############
# CONSTANTS #
#############
NA = "NA"
DEF_SIGN_CUTOFF = 0.05
DEF_THRESH = 25
DEF_EXON_LEN_NORM = 100.0

DEF_DPSI_THRESH = 5.0
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
    opt_parser.add_option("--in_prefix",
                          dest="in_prefix",
                          type="string",
                          help="""Prefix of output files created from
                                  createAS_CountTables. In createAS_CountTables,
                                  this is the -o option""",
                          default=None)
#   opt_parser.add_option("-i",
#                         dest="input_file",
#                         type="string",
#                         help="Resulting file from clusterASExons2.py",
#                         default=None)
#   opt_parser.add_option("--left_intron",
#                         dest="left_input",
#                         type="string",
#                         help="""Resulting file from clusterASExons2.py, which
#                                 contains the exclusion and inclusion counts
#                                 for just the left side of an intron retention
#                                 event.""",
#                         default=None)
#   opt_parser.add_option("--right_intron",
#                         dest="right_input",
#                         type="string",
#                         help="""Resulting file from clusterASExons2.py, which
#                                 contains the exclusion and inclusion counts
#                                 for just the right side of an intron retention
#                                 event.""",
#                         default=None)
#   opt_parser.add_option("--lenNormalized_counts",
#                         dest="lenNormalized_counts",
#                         type="string",
#                         help="""File containing length-normalized inclusion
#                                 exclusion counts. Used for PSI calculation,
#                                 not for statistcal significance.""",
#                         default=None)
#   opt_parser.add_option("--lenNormalized_left_intron",
#                         dest="lenNormalized_left_intron_counts",
#                         type="string",
#                         help="""File containing length-normalized
#                                 the left intron_retention counts.
#                                 Used for PSI calculation, not for
#                                 statistical significane.""",
#                         default=None)
#   opt_parser.add_option("--lenNormalized_right_intron",
#                         dest="lenNormalized_right_intron_counts",
#                         type="string",
#                         help="""File containing length-normalized
#                                 the right intron_retention counts.
#                                 Used for PSI calculation, not for
#                                 for statistical significance.""",
#                         default=None)
    opt_parser.add_option("--has_virtual",
                          dest="has_virtual",
                          action="store_true",
                          help="""Gives flags that a virtual reference is being
                                  used.""",
                          default=False)
    opt_parser.add_option("--jcn_seq_len",
                          dest="jcn_seq_len",
                          type="int",
                          help="""Junction length. Used as an option in
                                  getASEventReadCounts.py""",
                          default=None)
    opt_parser.add_option("--output_dir",
                          dest="output_dir",
                          type="string",
                          help="Directory to place output files.",
                          default=None)
    opt_parser.add_option("--out_prefix",
                          dest="prefix",
                          type="string",
                          help="Prefix of all output files. DEF=None",
                          default=None)
#   opt_parser.add_option("--psi_output_most_sign",
#                         dest="psi_output",
#                         type="string",
#                         help="""Output file that will contain the PSI values
#                                 for all events and samples that are
#                                 signficantly spliced.""",
#                         default=None)
#   opt_parser.add_option("--psi_output_sign_by_samp",
#                         dest="psi_output_by_samp",
#                         type="string",
#                         help="""Output file that will contain the PSI values
#                                 for all events and samples that are
#                                 signficantly differentially spliced where
#                                 multiple testing is not done for all samples
#                                 tested against the virtual reference""",
#                         default=None)
#   opt_parser.add_option("--all_psi_output",
#                         dest="all_psi_output",
#                         type="string",
#                         help="""Output file that will contain the PSI values
#                                 for all events and samples that pass minimum
#                                 count thresholds""",
#                         default=None)
#   opt_parser.add_option("--left_intron_all_psi_output",
#                         dest="left_intron_all_psi_output",
#                         type="string",
#                         help="""Output file that will contain the PSI values
#                                 for the left side of intron retention
#                                 samples. Not required, but used for dPSI
#                                 thresholds when taking all splice events.""",
#                         default=None)
#   opt_parser.add_option("--right_intron_all_psi_output",
#                         dest="right_intron_all_psi_output",
#                         type="string",
#                         help="""Output file that will contain the PSI values
#                                 for the right side of intron retention
#                                 samples. Not required, but used for dPSI
#                                 thresholds when taking all splice events.""",
#                         default=None)
#   opt_parser.add_option("--recalculate_ref_psi",
#                         dest="recalculate_ref_psi",
#                         action="store_true",
#                         help="""The reference PSI given in input tables
#                                 should be recalculated due to changes in
#                                 thresholding for minimum input between
#                                 length-normalized and raw counts.""",
#                         default=False)
#   opt_parser.add_option("--pval_output",
#                         dest="pval_output",
#                         type="string",
#                         help="""Output file that will associate the
#                                 unadjusted and adjusted p-values for all
#                                 pairs that were tested.""",
#                         default=None)
#   opt_parser.add_option("--event_sum",
#                         dest="event_sum",
#                         type="string",
#                         help="""Output file that will contain the sum of the
#                                 exclusion and inclusion counts for every
#                                 sample that was considered signifcantly
#                                 affected.""",
#                         default=None)
    opt_parser.add_option("--thresh",
                          dest="threshold",
                          type="int",
                          help="""Threshold for minimum number of total reads
                                  in an event. Default=%d""" % DEF_THRESH,
                          default=DEF_THRESH)
    opt_parser.add_option("--min_dpsi_threshold",
                          dest="dpsi_threshold",
                          type="float",
                          help="""Threshold for minimum delta PSI value between
                                  the sample with the smallest and largest PSI.
                                  Events with dPSI values below the threshold
                                  will not be tested or reported. Def=%.2f""" % DEF_DPSI_THRESH,
                          default=DEF_DPSI_THRESH)
    opt_parser.add_option("--method",
                          dest="method",
                          type="string",
                          help="""Correction Method: "BH" - Benjamini & Hochberg,
                                  "bonferroni".  Must select these strings as
                                  the option""",
                          default=None)
    opt_parser.add_option("--sign_cutoff",
                          dest="sign_cutoff",
                          type="float",
                          help="""Cutoff of corrected p-value significance.
                                  Default=%.2f""" % DEF_SIGN_CUTOFF,
                          default=DEF_SIGN_CUTOFF)
    opt_parser.add_option("--weights",
                          dest="weights",
                          type="string",
                          help="""Comma separated list of weights given in the
                                  order of the samples in the table. Weights are
                                  used to create a weighted median. Default is
                                  equal weight for all samples.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
#    opt_parser.check_required("-i")
#    opt_parser.check_required("--psi_output_most_sign")
#    opt_parser.check_required("--pval_output")
#    opt_parser.check_required("--event_sum")
    opt_parser.check_required("--method")
    opt_parser.check_required("--in_prefix")
    opt_parser.check_required("--out_prefix")
    opt_parser.check_required("--jcn_seq_len")

    in_prefix = options.in_prefix
    prefix = options.prefix

    try:
        input_file = open(in_prefix + "_AS_exclusion_inclusion_counts.txt")
    except:
        print """Cannot find expected file %s_AS_exclusion_inclusion_counts.txt.
                 Please check that the same options is given from
                 combine_createAS_CountTables""" % prefix
        opt_parser.print_help()
        sys.exit(1)

    left_input_file_name = in_prefix + "_left_intron_counts.txt"
    right_input_file_name = in_prefix + "_right_intron_counts.txt"
    sum_thresh = options.threshold

    sign_cutoff = options.sign_cutoff

    dpsi_thresh = options.dpsi_threshold

    left_input_file = None
    right_input_file = None
    if left_input_file_name is None:
        print "Warning: No intron retention file given as input.  Will not calculate IR events."
    else:
        left_input_file = open(left_input_file_name)
        right_input_file = open(right_input_file_name)

    output_dir = formatDir(options.output_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    out_prefix = "%s/%s" % (output_dir, options.prefix)
    
    psi_out = open("%s_most_sign_PSI.txt" % out_prefix, "w")
    pval_out = open("%s_pairs_p_val.txt" % out_prefix, "w")

    has_virtual = options.has_virtual

    # Optional output files
    psi_out_by_samp = open("%s_sign_by_samp_PSI.txt" % out_prefix, "w")
   
    all_psi_output = open("%s_allPSI.txt" % out_prefix, "w")

    left_all_psi_output = open("%s_left_intron_retention_allPSI.txt" % out_prefix, "w")

    right_all_psi_output = open("%s_right_intron_retention_allPSI.txt" % out_prefix, "w")

    jcn_seq_len = options.jcn_seq_len

    recalculate_ref_psi = False
    lenNormalized_counts_event2PSIs = None
    lenNormalized_counts_event2total_counts = None
#   if options.lenNormalized_counts:
#       if ((not options.lenNormalized_left_intron_counts) or 
#           (not options.lenNormalized_right_intron_counts)):
#           print "Need to specify all length-normalized count files."
#           opt_parser.print_help()
#           sys.exit(1)
#   
    recalculate_ref_psi = True
    lenNormalized_counts = open(in_prefix + "_AS_exclusion_inclusion_counts_lenNorm.txt")
    (lenNormalized_counts_event2total_counts,
     lenNormalized_counts_event2PSIs) = buildDicts(lenNormalized_counts)
    lenNormalized_counts.close()

    left_lenNormalized_counts_event2total_counts = None
    left_lenNormalized_counts_event2PSIs = None
#   if options.lenNormalized_left_intron_counts:
#       if ((not options.lenNormalized_counts) or 
#           (not options.lenNormalized_right_intron_counts)):
#           print "Need to specify all length-normalized count files."
#           opt_parser.print_help()
#           sys.exit(1)

    left_lenNormalized_counts = open(in_prefix + "_left_intron_counts_lenNorm.txt")
    (left_lenNormalized_counts_event2total_counts,
     left_lenNormalized_counts_event2PSIs) = buildDicts(left_lenNormalized_counts)
    left_lenNormalized_counts.close()

    right_lenNormalized_counts_event2total_counts = None
    right_lenNormalized_counts_event2PSIs = None
#   if options.lenNormalized_right_intron_counts:
#       if ((not options.lenNormalized_counts) or 
#           (not options.lenNormalized_left_intron_counts)):
#           print "Need to specify all length-normalized count files."
#           opt_parser.print_help()
#           sys.exit(1)

    right_lenNormalized_counts = open(in_prefix + "_right_intron_counts_lenNorm.txt")
    (right_lenNormalized_counts_event2total_counts,       
     right_lenNormalized_counts_event2PSIs) = buildDicts(right_lenNormalized_counts)
    right_lenNormalized_counts.close()

#    if options.lenNormalized_counts:
#       if not jcn_seq_len:
#           print "If length normalized counts are specified, need to give jcn_seq_len"
#           opt_parser.print_help()
#           sys.exit(1)

    weights = None
    if options.weights:
        weights = map(float, options.weights.split(","))

        # Use R limma package
        try:
            r.library("limma")
        except:
            print """In order to use weighted median, please install the limma package from Bioconductor: 
                     http://www.bioconductor.org/packages/release/bioc/html/limma.html"""
            print """In R:\nsource("http://bioconductor.org/biocLite.R")\nbiocLite("limma")"""

    event_sum = open("%s_event_sum.txt" % out_prefix, "w")

    method = options.method
    if method != "BH" and method != "bonferroni":
        print "Wrong method indicated."
        opt_parser.print_help()
        sys.exit(1)

    # {event_type:[pval]}
    event_type2pvals = {}

    # {event:(col1, col2):pval_idx}
    event2pairs2idx = {} 

    # Additional pval holders tested by each sample against the reference
    # {event_type:col:[pval]}
    event_type2col2pvals = {}

    # {event:col:pval_idx}
    event2col2idx = {}

    # {event:{col:psi}}
    event2col2psi = {}

    # {event:{col:sum_counts}}
    event2col2sum = {}

    # For weighted median
    col2weights = None

    header = None
    total_samples = None
    for line in input_file:
        line = formatLine(line)

        if line.startswith("#"):
            header = line
            line_list = line.split("\t")
            samples = line_list[11:]
            total_samples = len(samples)
            if weights:
                if len(weights) != total_samples-1:
                    print "Weights for every sample needs to be given"
                    opt_parser.print_help()
                    sys.exit(1)

                col2weights = {}
                for i in range(1,total_samples):
                    col2weights[i-1] = weights[i-1]
            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:11])
        counts = line_list[11:]

        # If the reference is NA, then do not calculate anything
        if counts[0] == NA:
            continue

        if has_virtual:
            # Cannot do a comparison when virtual reference is low expressed
            if lenNormalized_counts_event2total_counts[event][0] == NA:
                continue
        
        lenNormalized_psis = [None for i in range(len(counts))]
        if lenNormalized_counts_event2PSIs:
            try:
                lenNormalized_psis = lenNormalized_counts_event2PSIs[event]
            except:
                print "Warning: Can't find event in lenNormalized psis: %s" % event
                continue

        event_type = getEventType(event)
        if event_type not in event_type2pvals:
            event_type2pvals[event_type] = []
        if event_type not in event_type2col2pvals:
            event_type2col2pvals[event_type] = {}

        # Fill PSI dict
        for i in range(total_samples):
            (psi, sum_ct) = getPSI_sample_sum(counts[i], sum_thresh,
                                              lenNormalized_psis[i])
            if event in event2col2psi:
                event2col2psi[event][i] = psi
                event2col2sum[event][i] = sum_ct
            else:
                event2col2psi[event] = {i:psi}
                event2col2sum[event] = {i:sum_ct}

        # Only psis in event2col2psi that passed the sum_thresh will be
        # present, for ref psi will be calculated from the median of the
        # existing values
        if recalculate_ref_psi and has_virtual:
            adj_psi, adj_totalCount = recalculateRefPSI(event2col2psi[event],
                                                        lenNormalized_counts_event2total_counts[event],
                                                        col2weights)
            event2col2psi[event][0] = adj_psi
            lenNormalized_counts_event2total_counts[event][0] = adj_totalCount

        if dPSI(event2col2psi[event]) < dpsi_thresh:
            for j in range(1,total_samples):
                if event in event2pairs2idx:
                    event2pairs2idx[event][(0,j)] = NA
                else:
                    event2pairs2idx[event] = {(0,j):NA}
                if event in event2col2idx:
                    event2col2idx[event][j] = NA
                else:
                    event2col2idx[event] = {j:NA}

            continue

        # Calculate p-val for intron retention later
        if event_type == "intron_retention":
            continue

        # Do pairwise comparisons with first column
        [col1_excl, col1_incl] = map(int,counts[0].split(";"))
        if recalculate_ref_psi and has_virtual:
            # Need to also adjust relative counts based on new PSI
            col1_excl, col1_incl = adjustRefCounts(event,
                                                   jcn_seq_len,
                                                   lenNormalized_counts_event2total_counts[event][0],
                                                   float(event2col2psi[event][0]),
                                                   col1_excl,
                                                   col1_incl)

        for j in range(1,total_samples):

            if j not in event_type2col2pvals[event_type]:
                event_type2col2pvals[event_type][j] = []

            [col2_excl, col2_incl] = map(int,counts[j].split(";"))

            # Both samples have to be non-zero
            if belowThreshold(sum_thresh, col1_excl, col1_incl,
                              col2_excl, col2_incl):
                if event in event2pairs2idx:
                    event2pairs2idx[event][(0,j)] = NA
                else:
                    event2pairs2idx[event] = {(0,j):NA}
                if event in event2col2idx:
                    event2col2idx[event][j] = NA
                else:
                    event2col2idx[event] = {j:NA}

                continue

            cur_len = len(event_type2pvals[event_type])
            cur_len2 = len(event_type2col2pvals[event_type][j])

            if event in event2pairs2idx:
                event2pairs2idx[event][(0,j)] = cur_len
            else:
                event2pairs2idx[event] = {(0,j):cur_len}	

            if event in event2col2idx:
                event2col2idx[event][j] = cur_len2
            else:
                event2col2idx[event] = {j:cur_len2}

            raw_pval = robjects.r['fisher.test'](r.matrix(robjects.IntVector([col1_excl,
                                                                              col1_incl,
                                                                              col2_excl,
                                                                              col2_incl]),
                                                                              nrow=2))[0][0]

            event_type2pvals[event_type].append(raw_pval)

            updateDictOfLists(event_type2col2pvals[event_type], j, raw_pval)

    # Now calculate intron retention
    if left_input_file:
        left_events2counts = getIntronLeftRightCounts(left_input_file)
        right_events2counts = getIntronLeftRightCounts(right_input_file)
    else:
        left_events2counts = {}
        right_events2counts = {}

    if left_all_psi_output:
        left_all_psi_output.write(header + "\n")
    if right_all_psi_output:
        right_all_psi_output.write(header + "\n")

    for event in left_events2counts:
        if event not in right_events2counts:
            continue

        allPSI_elems_left = []
        allPSI_elems_right = []

        left_length = len(left_events2counts[event])
        right_length = len(right_events2counts[event])

        lenNormalized_left_psis = [None for i in range(left_length)]
        lenNormalized_right_psis = [None for i in range(right_length)]

        if left_lenNormalized_counts_event2PSIs:
            try:
                lenNormalized_left_psis = left_lenNormalized_counts_event2PSIs[event]
            except:
                print "Warning: Could not find event in left_lenNormalized psis: %s" % event
                continue
        if right_lenNormalized_counts_event2PSIs:
            try:
                lenNormalized_right_psis = right_lenNormalized_counts_event2PSIs[event]
            except:
                print "Warning: Could not find event in right_lenNormalized psis: %s" % event
                continue

        # Fill PSI dict
        for i in range(left_length):
            (psi, sum_ct) = getPSI_sample_sum(left_events2counts[event][i], sum_thresh,
                                              lenNormalized_left_psis[i])
            allPSI_elems_left.append(psi)

            try:
                (psi, sum_ct) = getPSI_sample_sum(right_events2counts[event][i], sum_thresh,
                                                  lenNormalized_right_psis[i])
            except:
                pdb.set_trace()
            allPSI_elems_right.append(psi)


#           # Adding left and right PSI values
#           if left_col2_excl + left_col2_incl < sum_thresh:
#               allPSI_elems_left.append(NA)
#           else:
#               allPSI_elems_left.append(getPSI(left_col2_excl, left_col2_incl,
#                                               lenNormalized_left_psis[j]))

#           if right_col2_excl + right_col2_incl < sum_thresh:
#               allPSI_elems_right.append(NA)
#           else:
#               allPSI_elems_right.append(getPSI(right_col2_excl,
#                                                right_col2_incl,
#                                                lenNormalized_right_psis[j]))

        # Only psis in event2col2psi that passed the sum_thresh will be
        # present, for ref psi will be calculated from the median of the
        # existing values
        if recalculate_ref_psi and has_virtual:
            allPSI_elems_left[0] = recalculateRefPSI_list(allPSI_elems_left,
                                                          col2weights)
            allPSI_elems_right[0] = recalculateRefPSI_list(allPSI_elems_right,
                                                           col2weights)

        if dPSI(allPSI_elems_left) < dpsi_thresh or dPSI(allPSI_elems_right) < dpsi_thresh:

            for j in range(1,left_length):
                if event in event2pairs2idx:
                    event2pairs2idx[event][(0,j)] = NA
                else:
                    event2pairs2idx[event] = {(0,j):NA}
                if event in event2col2idx:
                    event2col2idx[event][j] = NA
                else:
                    event2col2idx[event] = {j:NA}

            continue

        [left_col1_excl, left_col1_incl] = map(int,left_events2counts[event][0].split(";"))
        [right_col1_excl, right_col1_incl] = map(int,right_events2counts[event][0].split(";"))


        if left_col1_excl + left_col1_incl < sum_thresh:
            continue # the reference must have a PSI

        if right_col1_excl + right_col1_incl < sum_thresh:
            continue # the reference must have a PSI

        # Adjust ref counts based on PSI
        if recalculate_ref_psi and has_virtual:
            left_col1_excl, left_col1_incl = adjustRefCounts(event,
                                                             jcn_seq_len,
                                                             left_lenNormalized_counts_event2total_counts[event][0],
                                                             float(allPSI_elems_left[0]), 
                                                             left_col1_excl, 
                                                             left_col1_incl)

            right_col1_excl, right_col1_incl = adjustRefCounts(event,
                                                               jcn_seq_len,
                                                               right_lenNormalized_counts_event2total_counts[event][0],
                                                               float(allPSI_elems_right[0]), 
                                                               right_col1_excl, 
                                                               right_col1_incl)


        for j in range(1,total_samples):

            [left_col2_excl, left_col2_incl] = map(int,left_events2counts[event][j].split(";"))
            [right_col2_excl, right_col2_incl] = map(int,right_events2counts[event][j].split(";"))

            if j not in event_type2col2pvals["intron_retention"]:
                event_type2col2pvals["intron_retention"][j] = []

            # Both samples have to be non-zero
            if (belowThreshold(sum_thresh, left_col1_excl, left_col1_incl,
                              left_col2_excl, left_col2_incl) or
                belowThreshold(sum_thresh, right_col1_excl, right_col1_incl,
                               right_col2_excl, right_col2_incl)):
                if event in event2pairs2idx:
                    event2pairs2idx[event][(0,j)] = NA
                else:
                    event2pairs2idx[event] = {(0,j):NA}
                if event in event2col2idx:
                    event2col2idx[event][j] = NA
                else:
                    event2col2idx[event] = {j:NA}
                continue

            cur_len = len(event_type2pvals["intron_retention"])
            cur_len2 = len(event_type2col2pvals["intron_retention"][j])

            if event in event2pairs2idx:
                event2pairs2idx[event][(0,j)] = cur_len
            else:
                event2pairs2idx[event] = {(0,j):cur_len}	

            if event in event2col2idx:
                event2col2idx[event][j] = cur_len2
            else:
                event2col2idx[event] = {j:cur_len2}

            left_pval = robjects.r['fisher.test'](r.matrix(robjects.IntVector([left_col1_excl, 
                                                                               left_col1_incl,
                                                                               left_col2_excl,
                                                                               left_col2_incl]),
                                                           nrow=2))[0][0] 
            
            right_pval = robjects.r['fisher.test'](r.matrix(robjects.IntVector([right_col1_excl, 
                                                                                right_col1_incl,
                                                                                right_col2_excl,
                                                                                right_col2_incl]),
                                                           nrow=2))[0][0] 

            combined_pval = (left_pval + right_pval) - left_pval * right_pval

            event_type2pvals["intron_retention"].append(combined_pval)

            updateDictOfLists(event_type2col2pvals["intron_retention"], j, combined_pval)


        # All samples have been processed, now print to allPSI
        if left_all_psi_output:
            left_all_psi_output.write(event + "\t" +
                                      "\t".join(allPSI_elems_left) + "\n")
        if right_all_psi_output:
            right_all_psi_output.write(event + "\t" +
                                      "\t".join(allPSI_elems_right) + "\n")

    if left_all_psi_output:
        left_all_psi_output.close()
    if right_all_psi_output:
        right_all_psi_output.close()

    # All pairs have been evaluated, so now do multiple testing correction on
    # everything
    event_type2adjusted_pvals = {}
    event_type2col2adjusted_pvals = {}

    for event_type in event_type2pvals:
        event_type2adjusted_pvals[event_type] = robjects.r['p.adjust'](robjects.FloatVector(event_type2pvals[event_type]),
                                                                       method) 
    
    for event_type in event_type2col2pvals:
        event_type2col2adjusted_pvals[event_type] = {}
        for col in event_type2col2pvals[event_type]:
            event_type2col2adjusted_pvals[event_type][col] = robjects.r['p.adjust'](robjects.FloatVector(event_type2col2pvals[event_type][col]),
                                                                                    method)

    # Now go through all events and only consider those that are signficant
    psi_out.write(header + "\n")
    if psi_out_by_samp:
        psi_out_by_samp.write(header + "\n")
    if all_psi_output:
        all_psi_output.write(header + "\n")

    for event in event2pairs2idx:
        sign_cols = set([])
        sign_cols2 = set([])
        event_type = getEventType(event)

        for pair in event2pairs2idx[event]:
            this_idx = event2pairs2idx[event][pair]
            this_idx2 = event2col2idx[event][pair[1]]
            if this_idx == NA:
                continue

            outline = "%s\t%d\t%d\t%f" % (event,
                                          pair[0], pair[1],
                                          event_type2pvals[event_type][this_idx])
            if psi_out_by_samp:
                outline += "\t%f" % event_type2col2adjusted_pvals[event_type][pair[1]][this_idx2]

            outline += "\t%f\n" % event_type2adjusted_pvals[event_type][this_idx]
            pval_out.write(outline)

            if event_type2adjusted_pvals[event_type][this_idx] < sign_cutoff:
                sign_cols.add(pair[0])
                sign_cols.add(pair[1])

            if psi_out_by_samp:
                if event_type2col2adjusted_pvals[event_type][pair[1]][this_idx2] < sign_cutoff:
                    sign_cols2.add(pair[0])
                    sign_cols2.add(pair[1])

        # Write out PSI for any significant samples
        # Significant across all samples
        if sign_cols != set([]):        
            psi_vals = []
            for i in range(total_samples):
                if i in sign_cols:
                    psi_vals.append(event2col2psi[event][i])
                else:
                    psi_vals.append(NA) 

            outline = "%s\t%s\n" % (event, 
                                    "\t".join(psi_vals))
            psi_out.write(outline)
   
        # Significant by samples 
        if sign_cols2 != set([]):        
            psi_vals = []
            for i in range(total_samples):
                if i in sign_cols2:
                    psi_vals.append(event2col2psi[event][i])
                    if event_sum:
                        event_sum.write("%s\t%d\t%s\n" % (event,
                                                          i,
                                                          event2col2sum[event][i]))
                else:
                    psi_vals.append(NA) 

            outline = "%s\t%s\n" % (event, 
                                    "\t".join(psi_vals))
            psi_out_by_samp.write(outline)

        # Print all psi
        if all_psi_output:
            psi_vals = []
            for i in range(total_samples):
                try:
                    psi_vals.append(event2col2psi[event][i])
                except:
                    psi_vals.append(NA)

            outline = "%s\t%s\n" % (event, 
                                    "\t".join(psi_vals))

            all_psi_output.write(outline)
    
    psi_out.close()
    psi_out_by_samp.close()
    all_psi_output.close()
    pval_out.close()

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def adjustRefCounts(event_str, jcn_seq_len, lengthNorm_total, ref_psi, excl, incl):
    raw_total = excl + incl
    ref_psi_proportion = ref_psi/100

    inclusion_isoform_len = None
    if "jcn_only" in event_str:
        # These events were not length normalized since each isoform has the same length
        inclusion_isoform_len = DEF_EXON_LEN_NORM 
    elif "alternative_donor" in event_str or "alternative_acceptor" in event_str:
        # This can be removed once bug in reporting inclusion regions are
        # fixed, then getInclIsoformLen can simply be used
        inclusion_isoform_len = getAA_ADInclIsoformLen(event_str, jcn_seq_len)
    else:
        inclusion_isoform_len = getInclIsoformLen(event_str, jcn_seq_len)

    if not inclusion_isoform_len:
        print "Error in obtaining isoform length from: %s" % event_str
        sys.exit(1)

    adj_incl = int(round((float(inclusion_isoform_len)/DEF_EXON_LEN_NORM) *
                          ref_psi_proportion * lengthNorm_total))

    # Do to rounding, do not want to have negative values
    if adj_incl > raw_total:
        adj_incl = raw_total

    adj_excl = raw_total - adj_incl

    return adj_excl, adj_incl

def belowThreshold(sum_thresh, col1_excl, col1_incl, col2_excl, col2_incl):

    if col1_excl + col1_incl < sum_thresh:
        return True

    if col2_excl + col2_incl < sum_thresh:
        return True

    return False

def buildDicts(lenNormalized_counts_file):
    lenNormalized_counts_event2total_counts = {}
    lenNormalized_counts_event2psis = {}
    for line in lenNormalized_counts_file:
        line = formatLine(line)
        if line.startswith("#"):
            continue

        line_list = line.split("\t")

        event = "\t".join(line_list[0:11])
        counts_list = line_list[11:]
        
        lenNormalized_counts_event2total_counts[event] = []
        lenNormalized_counts_event2psis[event] = []

        for counts in counts_list:
            if counts == NA:
                lenNormalized_counts_event2psis[event].append(NA)
                lenNormalized_counts_event2total_counts[event].append(NA)
                continue

            [excl, incl] = map(int,counts.split(";"))
            total = excl + incl 
            if total == 0:
                lenNormalized_counts_event2psis[event].append(NA)
            else:
                psi_val = (float(incl)/total) * 100
                lenNormalized_counts_event2psis[event].append("%.2f" % psi_val)

            lenNormalized_counts_event2total_counts[event].append(total)
   
    return lenNormalized_counts_event2total_counts, lenNormalized_counts_event2psis 
    
def dPSI(psi_vals):
    minPSI = 130.0
    maxPSI = -10.0

    if type(psi_vals) == type({}):
        this_psi_vals = psi_vals.values()
    else:
        this_psi_vals = list(psi_vals)

    for psi_str in this_psi_vals:
        if psi_str == NA:
            continue
        psi = float(psi_str)

        if psi < minPSI:
            minPSI = psi
        if psi > maxPSI:
            maxPSI = psi

    return abs(maxPSI - minPSI)
    

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getAA_ADInclIsoformLen(event_str, jcn_seq_len):
    event_str_list = event_str.split("\t")

    excl_jcns = event_str_list[5].split(";")
    
    chr, incl_start, incl_end = convertCoordStr(event_str_list[6])

    alt_start_or_end = determineAltStartOrEnd(incl_start, incl_end,
                                                  excl_jcns)

    ordered_pos, inclusion_pos_idx = getSSOrder(alt_start_or_end, incl_start, incl_end, excl_jcns)
     
    isoform_lengths = getAD_AA_isoform_lengths(alt_start_or_end,
                                               ordered_pos, jcn_seq_len)

    return isoform_lengths[inclusion_pos_idx]

def getInclIsoformLen(event_str, jcn_seq_len):
    """
    Length of the inclusion isoform equals sum of junction lengths + exon
    lengths
    """
    event_str_list = event_str.split("\t")
    jcns = elems_split(event_str_list[6])
    jcns.extend(elems_split(event_str_list[9]))

    exon_len = 0
    exons = elems_split(event_str_list[8])
    for exon in exons:
        if exon == "" or exon == "None":
            continue
        chr, exon_start, exon_end = convertCoordStr(exon)        

        this_len = exon_end - exon_start + 1 
        exon_len += this_len

    return exon_len + (jcn_seq_len * (len(jcns)))
    

def getPSI_sample_sum(excl_incl_ct_str, sum_thresh, lenNormalized_PSI):
    if excl_incl_ct_str == NA:
        return NA, NA

    excl_str, incl_str = excl_incl_ct_str.split(";")

    excl = float(excl_str)
    incl = float(incl_str)

    if (excl + incl) <  sum_thresh:  
        return (NA, NA)
   
    total_ct = "%d" % int(excl + incl)

    if lenNormalized_PSI:
        return lenNormalized_PSI, total_ct

    psi = (incl/(incl + excl)) * 100

    psi_str = "%.2f" % psi

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

def getPSI(excl, incl, lenNormalizedPSI):
    if lenNormalizedPSI:
        return lenNormalizedPSI
    excl = float(excl)
    incl = float(incl)
    return "%.2f" % ((incl/(incl + excl)) * 100)
    
def getSSOrder(alt_start_or_end, inclusion_start, inclusion_end,
               exclusion_str_list):

    unordered_pos = []
    for_ordered_inclusion_pos = None

    if alt_start_or_end == "alt_start":
        unordered_pos.append(inclusion_start)
        for_ordered_inclusion_pos = inclusion_start
        for exclusion_str in exclusion_str_list:
            chr, start, end = convertCoordStr(exclusion_str)
            unordered_pos.append(start)
    else:
        unordered_pos.append(inclusion_end)
        for_ordered_inclusion_pos = inclusion_end
        for exclusion_str in exclusion_str_list:
            chr, start, end = convertCoordStr(exclusion_str)
            unordered_pos.append(end)

    ordered_pos = list(unordered_pos)
    ordered_pos.sort()

    return ordered_pos, ordered_pos.index(for_ordered_inclusion_pos)

def elems_split(jcn_str):
    jcn_list1 = jcn_str.split(";")
    returnJcns = []
    for jcn_elem in jcn_list1:
        if jcn_elem == "":
            continue
        jcn_list2 = jcn_elem.split(",") 
        for jcn in jcn_list2:
            if jcn == "":
                continue
            returnJcns.append(jcn)

    return returnJcns
    
def recalculateRefPSI(col2psi_str, col2total_count, col2weights):
    """
    Reference is assumed to be col 0
    """
    psi_vals = []
    total_vals = []
    weights = []
    for col in col2psi_str:
        if col == 0:
            continue
        if col2psi_str[col] == NA:
            continue
        psi_vals.append(float(col2psi_str[col]))
        total_vals.append(col2total_count[col])
        if col2weights:
            weights.append(col2weights[col-1])

    if col2weights:
        median_psi = r['weighted.median'](robjects.FloatVector(psi_vals), 
                                          robjects.FloatVector(weights))[0]
        median_total = int(round(r['weighted.median'](robjects.IntVector(total_vals),                                
                                                      robjects.FloatVector(weights))[0]))
    else:
        median_psi = r['median'](robjects.FloatVector(psi_vals))[0]
        try:
            median_total = int(round(r['median'](robjects.IntVector(total_vals))[0]))
        except:
            pdb.set_trace()

    return "%.2f" % median_psi, median_total


def recalculateRefPSI_list(psi_list, col2weights):
    """
    Reference is assumed to be col 0
    """
    vals = []
    weights = []
    for i in range(1,len(psi_list)):
        if psi_list[i] == NA:
            continue
        vals.append(float(psi_list[i]))
        if col2weights:
            weights.append(col2weights[i-1])

    if col2weights:
        median_psi = r['weighted.median'](robjects.FloatVecotr(vals),
                                          robjects.FloatVector(weights))[0]
    else:
        median_psi = r['median'](robjects.FloatVector(vals))[0]

    return "%.2f" % median_psi
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
