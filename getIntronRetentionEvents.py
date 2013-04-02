#!/lab/64/bin/python
# getIntronRetentionEvents.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Takes files generated from getASEventReadCounts and fisherTest_exclusionInclusion.py
   for the intron retention events.
"""

import sys
import optparse 
import math
import pdb
import os
import pickle

import rpy2.robjects as robjects


from getASEventReadCounts import isAnnotated, getAnnotatedIntronCoords, getAnnotatedExonCoords, getAllEventStr, findAdjacentSharedRegion
r = robjects.r
#############
# CONSTANTS #
#############
DEF_PVAL_CUTOFF = 0.05
BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
COORD_COUNT = "%s/coordReadCounts.py" % BIN_DIR
if not os.path.exists(COORD_COUNT):
    print "ERROR: coordReadCounts.py needs to be in the same directory."
    sys.exit(1)

HOST = "localhost"
USER = "root"
PASSWD = ""
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
    opt_parser.add_option("-l",
                          dest="left_file",
                          type="string",
                          help="""File from left side of intron.  Must include
                                  p-vals.""",
                          default=None)
    opt_parser.add_option("-r",
                          dest="right_file",
                          type="string",
                          help="""File from right side of intron.  Must include
                                  p-vals.""",
                          default=None)
    opt_parser.add_option("-p",
                          dest="prefix",
                          type="string",
                          help="""File prefix for output files.  Should be the
                                  same as the one used for
                                  getASEventReadCounts.py""",
                          default=None)
#   opt_parser.add_option("--sqlite_db_dir",
#                         dest="sqlite_db_dir",
#                         type="string",
#                         help="""Location of sqlite database. If sqlite database
#                                 is used, will override usage of MySQL database.""",
#                         default=None)
#   opt_parser.add_option("--host",
#                         dest="host",
#                         type="string",
#                         help="MySQL database host. Def=\'%s\'" % HOST,
#                         default=HOST)
#   opt_parser.add_option("--user",
#                         dest="user",
#                         type="string",
#                         help="MySQL database user. Def=\'%s\'" % USER,
#                         default=USER)
#   opt_parser.add_option("--passwd",
#                         dest="passwd",
#                         type="string",
#                         help="MySQL database password. Def=\'%s\'" % PASSWD,
#                         default=PASSWD)
#   opt_parser.add_option("-d",
#                         dest="as_db",
#                         type="string",
#                         help="Database of annotated AS events.",
#                         default=None)
    opt_parser.add_option("--intron",
                          dest="intron_pk_file",
                          type="string",
                          help="Pickle file containing intron coordinates.",
                          default=None)
    opt_parser.add_option("--exon",
                          dest="exon_pk_file",
                          type="string",
                          help="Pickle file containing exon coordinates.",
                          default=None)
    opt_parser.add_option("--all_as_event",
                          dest="all_as_event_file",
                          type="string",
                          help="""File containing all information for all AS
                                  events.""",
                          default=None)
    opt_parser.add_option("--lengthNorm",
                          dest="lengthNorm",
                          action="store_true",
                          help="""Indicates that read counts in the
                                  all_as_event file should be normalized by the length""",
                          default=False)
#   opt_parser.add_option("--method",
#                         dest="method",
#                         type="string",
#                         help="""Type of correction method:
#                                 'BH' - Benjamini & Hochberg,
#                                 'bonferroni'""",
#                         default=None)
    opt_parser.add_option("--by_chr",
                          dest="this_chr",
                          type="string",
                          help="""If running samples by each chromosome, the
                                  chromosome is given in this option.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-l")
    opt_parser.check_required("-r")
#    opt_parser.check_required("-d")
#    opt_parser.check_required("--method")
    opt_parser.check_required("--all_as_event")
    opt_parser.check_required("--intron")
    opt_parser.check_required("--exon")


#    method = options.method

#   if method != "BH" and method != "bonferroni":
#       print "Wrong method given."
#       opt_parser.print_help()
#       sys.exit(1)
    
#    as_db = options.as_db
    lengthNorm = options.lengthNorm

    this_chr = options.this_chr
    
#   db = None
#   if options.sqlite_db_dir:
#       import pysqlite_wrap
#       db = pysqlite_wrap.DB(options.sqlite_db_dir)
#   else: # Use MySQL database
#       import mysqldb_wrap
#       db = mysqldb_wrap.DB(options.host, options.user, options.passwd)

    left_file = open(options.left_file)
    right_file = open(options.right_file)

    all_as_event_file = open(options.all_as_event_file, "w")

    prefix = options.prefix
    if prefix:
        # Create output files
        left_no_partner_out = open(prefix + "_IR_left_no_partner.txt", "w")
        right_no_partner_out = open(prefix + "_IR_right_no_partner.txt", "w")
        diff_direction_out = open(prefix + "_IR_diff_direction.txt", "w")
#        affected_out = open(prefix + "_IR_events.txt", "w")    
    else:
        left_no_partner_out = open("IR_left_no_partner.txt", "w")
        right_no_partner_out = open("IR_right_no_partner.txt", "w")
        diff_direction_out = open("IR_diff_direction.txt", "w")
#        affected_out = open("IR_events.txt", "w")    

    # {chr: set[(start, end, strand)])
#    annotated_introns = getAnnotatedIntronCoords(db, as_db, this_chr)
    intron_pk_file = open(options.intron_pk_file)
    annotated_introns = pickle.load(intron_pk_file)
    intron_pk_file.close()

#   (annotated_exons,
#    annotated_internal_exons,
#    annotated_exons_no_strand,
#    annotated_exons_by_strand,
#    annotated_exon_search_tree) = getAnnotatedExonCoords(db, as_db, this_chr)
    exon_pk_file = open(options.exon_pk_file)
    annotated_exons_by_strand = pickle.load(exon_pk_file)
    exon_pk_file.close()

    left_line_list = left_file.readlines()
    right_line_list = right_file.readlines()

    left_file.close()
    right_file.close()

    left_intron2line = {}
    right_intron2line = {}
    
    # Building Left Dict
    for i in range(len(left_line_list)):
        introns = parseLine(left_line_list[i])
    
        for intron in introns:
            left_intron2line[intron] = i

    # Building Right Dict
    for i in range(len(right_line_list)):
        introns = parseLine(right_line_list[i])

        for intron in introns:
            right_intron2line[intron] = i


    combined_pvals = []
    pval_intron_tuples = []

    # Processing, starting from left dict
    for intron in left_intron2line:
        if intron not in right_intron2line:
            left_no_partner_out.write(left_line_list[left_intron2line[intron]])
        elif isInDiffDirection(intron, left_intron2line, right_intron2line,
                               left_line_list, right_line_list):
            combined_line = combineLines(left_line_list[left_intron2line[intron]],
                                         right_line_list[right_intron2line[intron]])

            diff_direction_out.write(combined_line + "\t" + intron + "\n")

            # Print all events to allEvent file
            out_str = getAllEventInfoLine(combined_line, intron, annotated_exons_by_strand, annotated_introns, lengthNorm)
            all_as_event_file.write(out_str + "\n")
        else:
#           combined_pval = getCombinedPval(left_line_list[left_intron2line[intron]],
#                                           right_line_list[right_intron2line[intron]])

            # Write all IR events to allEventInfo file.
            combined_line = combineLines(left_line_list[left_intron2line[intron]],
                                         right_line_list[right_intron2line[intron]])

            # Only when both samples have counts, will a p_val be associated.
#           if isValidIR(combined_line):
#               combined_pvals.append(combined_pval)

#               pval_intron_tuples.append((combined_pval, intron))

            # Print all events to allEvent file
            out_str = getAllEventInfoLine(combined_line, intron, annotated_exons_by_strand, annotated_introns, lengthNorm)
            all_as_event_file.write(out_str + "\n")

    # Checking for right_no_partner
    for intron in right_intron2line:
        if intron not in left_intron2line:
            right_no_partner_out.write(right_line_list[right_intron2line[intron]])
            

# No longer doing p-value test here
#   # Adjust the pvalues
#   adj_pvals_rVec = robjects.r['p.adjust'](robjects.FloatVector(combined_pvals), 
#                      method)

#   adj_pvals = []

#   for p_val in adj_pvals_rVec:
#       adj_pvals.append(p_val)

#   adj_pvals_sort = list(adj_pvals)
#   adj_pvals_sort.sort()
#  
#   for pval in adj_pvals_sort:
#       for i in range(len(pval_intron_tuples)):
#           if pval_intron_tuples[i] == None:
#               continue
#           this_pval = adj_pvals[i]
#           intron = pval_intron_tuples[i][1]
#           if this_pval == pval:
#               combined_line = combineLines(left_line_list[left_intron2line[intron]],
#                                            right_line_list[right_intron2line[intron]])

#               affected_out.write(combined_line + "\t" + intron + "\t" +
#                                  repr(pval_intron_tuples[i][0]) + "\t" + 
#                                  repr(pval) + "\n")

#               pval_intron_tuples[i] = None
#               break
#           
#   affected_out.close()

    all_as_event_file.close()

#   # Add Constitutive counts
#   # Produce exon coord file
#   exon_coord_file = open("tmp2_exon_coord_file.txt", "w")

#   for coord in const_exons:

#       out_coord = coord.replace("_", "\t")
#       exon_coord_file.write(out_coord + "\n")

#       exon_coord_file.close()

#       mapped_file1_name = "tmp2_exon_coord_file1_readCounts.txt"
#       mapped_file2_name = "tmp2_exon_coord_file2_readCounts.txt"

#       cmd = "python %s -l %d --coords tmp2_exon_coord_file.txt " % (COORD_COUNT,
#                                                                    read_length)

#       # Run First Set of Reads through coordReadCounts
#       first_cmd = cmd + "--reads %s -o %s" % (genome_read_file1,
#                                               mapped_file1_name)

#       runCmd(first_cmd, SHELL, True)


    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def combineLines(left_line, right_line):
    """ 
    Returns a combination of the two lines
    """
    left_line_formatted = formatLine(left_line)
    right_line_formatted = formatLine(right_line)

    return left_line_formatted + "\t" + right_line_formatted

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getAllEventInfoLine(combined_line, intron, annotated_exons, annotated_introns, lengthNorm):
    (e_or_i, gene_name, chr, strand,
     left_intron_start, right_intron_end,
     excl_cts_raw, excl_cts_lenNorm,
     ie_cts_raw, ie_cts_lenNorm) = getEventInfo(combined_line)

    ie_jcns = ["%s:%d-%d" % (chr, left_intron_start - 1, left_intron_start)]
    ie_jcns.append("%s:%d-%d" % (chr, right_intron_end, right_intron_end + 1))
 
    const_regions = [] 
    const_str = findAdjacentSharedRegion(chr, strand, annotated_exons,
                                         left_intron_start - 1, "P")
    if const_str:
        const_regions.append(const_str)
    const_str = findAdjacentSharedRegion(chr, strand, annotated_exons,
                                         right_intron_end + 1, "N")
    if const_str:
        const_regions.append(const_str)

    ie_ct_list_raw = map(int, ie_cts_raw.split(";"))

    ie_ct_list_lenNorm = map(int, ie_cts_lenNorm.split(";"))

    ie_ct_raw = 0
    for ct in ie_ct_list_raw:
        ie_ct_raw += ct
    ie_ct_lenNorm = 0
    for ct in ie_ct_list_lenNorm:
        ie_ct_lenNorm += ct

    # Now print to all AS Event string
    label = "K"
    if not isAnnotated(intron, annotated_introns):
        label = "N"

    out_str = getAllEventStr(label,
                             "intron_retention",
                             gene_name,
                             chr,
                             strand,
                             intron,
                             "", # Inclusion junctions
                             "", # Exclusion exons
                             "", # Inclusion exons
                             ";".join(ie_jcns),
                             ";".join(const_regions),
                             excl_cts_raw,
                             "",
                             excl_cts_raw,
                             None,
                             "",
                             "",
                             None,
                             None,
                             ie_cts_raw,
                             ie_ct_raw,
                             "",
                             None)

    return out_str


def getEventInfo(combined_line):
    """
    Takes the combined information from both lines and output the information
    needed.
    e_or_i, gene_name, chr, strand,
    left_intron_start, right_intron_end
 
                (e_or_i, gene_name, chr, strand,
                 left_intron_start, right_intron_end,
                 excl_cts_raw, excl_cts_lenNorm,
                 ie_cts_raw, ie_cts_lenNorm) = getEventInfo(combined_line)

    """
    line_elems = combined_line.split("\t")

    e_or_i = line_elems[0]
    gene_name = line_elems[1]
    chr = line_elems[2]

    strand_first = line_elems[3]    
#    strand_second = line_elems[15]
    strand_second = line_elems[13]

    if strand_first != strand_second:
        strand = "."
    else:
        strand = strand_first

    left_intron_start = int(line_elems[4])
#    right_intron_end = int(line_elems[16])
    right_intron_end = int(line_elems[14])

    excl_cts_raw = int(line_elems[6])
    excl_cts_lenNorm = int(line_elems[8])

    ie_left_raw = None
    ie_right_raw = None
    ie_left_lenNorm = None
    ie_right_lenNorm = None

#    if lengthNorm: # Each ie counts need to be renormalized considering both ends
    ie_left_raw = int(line_elems[7])
#        ie_right_raw = int(round(float(line_elems[19])/2))
    ie_right_raw = int(line_elems[17])
    ie_left_lenNorm = int(round(float(line_elems[9])/2))
#        ie_right_lenNorm = int(round(float(line_elems[21])/2))
    ie_right_lenNorm = int(round(float(line_elems[19])/2))
#   else:
#       ie_left_raw = int(line_elems[7])
#         ie_right_raw = int(line_elems[19])
#       ie_right_raw = int(line_elems[17])
#       ie_left_lenNorm = int(line_elems[9])
#         ie_right_lenNorm = int(line_elems[21])
#       ie_right_lenNorm = int(line_elems[19])

    ie_cts_raw = "%d;%d" % (ie_left_raw,
                          ie_right_raw)
    ie_cts_lenNorm = "%d;%d" % (ie_left_lenNorm,
                          ie_right_lenNorm)

    return (e_or_i, gene_name, chr, strand, left_intron_start,
            right_intron_end, excl_cts_raw, excl_cts_lenNorm,
            ie_cts_raw, ie_cts_lenNorm)
        
def getCombinedPval(left_line, right_line):
    left_line = formatLine(left_line)
    right_line = formatLine(right_line)

    left_pval = float(left_line.split("\t")[-2])
    right_pval = float(right_line.split("\t")[-2])

    return (left_pval + right_pval) - left_pval * right_pval


def getEorI(lines, i):
    line = lines[i]

    line = formatLine(line)

    return line.split("\t")[0] 

def isInDiffDirection(intron, left_intron2line, right_intron2line,
                      left_line_list, right_line_list):
    left_e_or_i = getEorI(left_line_list, left_intron2line[intron])
    right_e_or_i = getEorI(right_line_list, right_intron2line[intron])

    if left_e_or_i != right_e_or_i:
        return True

    return False

def isValidIR(combined_line):
    """
    If either Sample 1 or Sample 2 both have zero counts, then it is not
    considered valid.
    """
    line_list = combined_line.split("\t")

    # Sample 1 exclusion counts
    raw_excl = int(line_list[6])
#    raw_excl += int(line_list[18])
    # No longer doing p-val calculation here.
    raw_excl += int(line_list[16])

    raw_incl = int(line_list[7])
#    raw_incl += int(line_list[19])
    raw_incl += int(line_list[17])

    lenNorm_excl = int(line_list[8])
#    lenNorm_excl += int(line_list[20])
    lenNorm_excl += int(line_list[18])
    
    lenNorm_incl = int(line_list[9])
#    lenNorm_incl += int(line_list[21])
    lenNorm_incl += int(line_list[19])

    if raw_excl == 0 and raw_incl == 0:
        return False
    if lenNorm_excl == 0 and lenNorm_incl == 0:
        return False

    return True
    

def parseLine(line):
    line = formatLine(line)

    line_list = line.split("\t")
    
    intron_elem = line_list[5]

    return intron_elem.split(",")
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
