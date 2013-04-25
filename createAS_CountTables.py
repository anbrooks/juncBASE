#!/lab/64/bin/python
# clusterASExons2.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Uses the *_all_AS_event_info.txt to create a full summary table. 
"""

import sys
import optparse 
import pdb
import os

from compareSampleSets import getSamples
from getASEventReadCounts import normalizeByLen

#############
# CONSTANTS #
#############
NA = "NA"

DEF_EXON_LEN_NORM = 100.0
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
                          dest="root_dir",
                          type="string",
                          help="""Root directory that contains subdirectories
                                  with output from getASEventReadCounts""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_file_prefix",
                          type="string",
                          help="""Output files that will contain all exclusion
                                  inclusion counts for every sample. As well as
                                  files for intron retention calculation.
                                  Finally, length-normalized counts are also
                                  produced.""",
                          default=None)
    opt_parser.add_option("--jcn_seq_len",
                          dest="jcn_seq_len",
                          type="int",
                          help="""Value used in getASEventReadCounts""", 
                          default=None)
#   opt_parser.add_option("--left_intron",
#                         dest="left_intron_file",
#                         type="string",
#                         help="""Output file that will contain the left side
#                                 of intron retention events.  Significant
#                                 intron retention events are identified in a
#                                 different way.""",
#                         default=None)
#   opt_parser.add_option("--right_intron",
#                         dest="right_intron_file",
#                         type="string",
#                         help="""Output file that will contain the right side
#                                 of intron retention events.  Significant
#                                 intron retention events are identified in a
#                                 different way.""",
#                         default=None)
#   opt_parser.add_option("--psi_output_file",
#                         dest="psi_output_file",
#                         type="string",
#                         help="""Optional: Output file that will contain the PSI value
#                                 for every event. It corresponds to the
#                                 output_file.""",
#                         default=None)
#   opt_parser.add_option("--psi_left_intron_file",
#                         dest="psi_left_intron_file",
#                         type="string",
#                         help="""Optional: Output file that will contain the PSI value
#                                 for the left side of intron retention
#                                 events.""",
#                         default=None)
#   opt_parser.add_option("--psi_right_intron_file",
#                         dest="psi_right_intron_file",
#                         type="string",
#                         help="""Optional: Output file that will contain the PSI value
#                                 for the right side of intron retention
#                                 events.""",
#                         default=None)
    opt_parser.add_option("-s",
                          dest="samples",
                          type="string",
                          help="""Comma separated list of the samples that will
                                  be used or a file of sample names.  The order which they are given is
                                  the order in the output of the file.""",
                          default=None)
#   opt_parser.add_option("--lengthNorm",
#                         dest="lengthNorm",
#                         action="store_true",
#                         help="""Flag to indicate length normalization was
#                                 done on the counts. Used for splitting the IR
#                                 counts back into left and right counts""",
#                         default=False)
    opt_parser.add_option("--which_chr",
                          dest="which_chr",
                          type="string",
                          help="""When running by chromsome, it will find the
                                  appropriate files given the expected directory
                                  structure.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-d")
    opt_parser.check_required("-o")
    opt_parser.check_required("-s")
    opt_parser.check_required("--jcn_seq_len")
#   opt_parser.check_required("--left_intron")
#   opt_parser.check_required("--right_intron")

    root_dir = options.root_dir

    if os.path.exists(root_dir):		
        root_dir = os.path.abspath(root_dir)
    else:
        print "Root directory does not exist: %s" % root_dir
        opt_parser.print_help()
        sys.exit(1)

    if not root_dir.endswith("/"):
        root_dir += "/"
    
    prefix = options.output_file_prefix
    
    output_file_name = prefix + "_AS_exclusion_inclusion_counts.txt"
    left_file_name = prefix + "_left_intron_counts.txt"
    right_file_name = prefix + "_right_intron_counts.txt"

    lenNorm_output_file_name = prefix + "_AS_exclusion_inclusion_counts_lenNorm.txt"
    lenNorm_left_file_name = prefix + "_left_intron_counts_lenNorm.txt"
    lenNorm_right_file_name = prefix + "_right_intron_counts_lenNorm.txt"

    samples = getSamples(options.samples)
    
    num_samples = len(samples)

    which_chr = None
    if options.which_chr:
        which_chr = options.which_chr

    jcn_seq_len = options.jcn_seq_len

    # Will equal none if no length normalization occurred
#    lengthNorm = options.lengthNorm

    # {event:sample:(excl, incl)}
    event2sample2raw_counts = {}
    event2sample2lenNorm_counts = {}

    # {event:(set(genes), strand)
    event2genesStrand = {}    

    # {event:sample:(excl, incl)}
    left_intron2sample2raw_counts = {}
    right_intron2sample2raw_counts = {}
    left_intron2sample2lenNorm_counts = {}
    right_intron2sample2lenNorm_counts = {}

    first_sample = samples.pop(0)

    # Initialize dictionary with first sample
    if which_chr:
        first_file = open(root_dir + first_sample + "/" + first_sample + "_" + which_chr + 
                          "/" + first_sample + "_" + which_chr + "_all_AS_event_info.txt")
    else:
        first_file = open(root_dir + first_sample + "/" + first_sample + "_all_AS_event_info.txt")

    for line in first_file:
        line = formatLine(line)

        event_key, raw_count_str, lenNorm_count_str, genes, strand, const_region = getKeyandCount(line)

        event2sample2raw_counts[event_key] = {first_sample:raw_count_str}
        event2sample2lenNorm_counts[event_key] = {first_sample:lenNorm_count_str}
        
        updateGeneStrand(event2genesStrand, event_key, genes, strand, const_region)

        if "intron_retention" in line:
            left_raw_count_str, right_raw_count_str = getIRKeyandCounts(line, jcn_seq_len, False)
            left_lenNorm_count_str, right_lenNorm_count_str = getIRKeyandCounts(line, jcn_seq_len, True)
            left_intron2sample2raw_counts[event_key] = {first_sample:left_raw_count_str}
            right_intron2sample2raw_counts[event_key] = {first_sample:right_raw_count_str}
            left_intron2sample2lenNorm_counts[event_key] = {first_sample:left_lenNorm_count_str}
            right_intron2sample2lenNorm_counts[event_key] = {first_sample:right_lenNorm_count_str}

    first_file.close()
           
    # Now populate with the rest of the files.
    for samp in samples:
        if which_chr:
            samp_file = open(root_dir + samp + "/" + samp + "_" + which_chr +
                             "/" + samp + "_" + which_chr + "_all_AS_event_info.txt") 
        else:
            samp_file = open(root_dir + samp + "/" + samp + "_all_AS_event_info.txt") 
        for line in samp_file:
            line = formatLine(line)

            event_key, raw_count_str, lenNorm_count_str, genes, strand, const_region = getKeyandCount(line)

            if event_key not in event2sample2raw_counts:
                print "Event from sample, not in dict: %s, %s" % (samp,
                                                                  event_key)
                continue

            event2sample2raw_counts[event_key][samp] = raw_count_str
            event2sample2lenNorm_counts[event_key][samp] = lenNorm_count_str

            updateGeneStrand(event2genesStrand, event_key, genes, strand, const_region)
        
            if "intron_retention" in line:
                left_raw_count_str, right_raw_count_str = getIRKeyandCounts(line, jcn_seq_len, False)
                left_lenNorm_count_str, right_lenNorm_count_str = getIRKeyandCounts(line, jcn_seq_len, True)
                left_intron2sample2raw_counts[event_key][samp] = left_raw_count_str
                right_intron2sample2raw_counts[event_key][samp] = right_raw_count_str
                left_intron2sample2lenNorm_counts[event_key][samp] = left_lenNorm_count_str
                right_intron2sample2lenNorm_counts[event_key][samp] = right_lenNorm_count_str

        samp_file.close()

    # Now print out
    output_file = open(output_file_name, "w")
    lenNorm_file = open(lenNorm_output_file_name, "w")

#   psi_output_file = None
#   if psi_output_file_name:
#       psi_output_file = open(psi_output_file_name, "w")

    skipped_ir_events = set([])

    # Print header
    header = "#Contains_Novel_or_Only_Known(Annotated)_Junctions\tas_event_type\tgene_name\tchr\tstrand\t"
    header += "exclusion_junctions\tinclusion_junctions\texclusion_exons\t"
    header += "inclusion_exons\tintron-exon_junctions\tneighboring_constitutive_exons\t"
    header += "%s\t%s\n" % (first_sample,
                            "\t".join(samples))
    output_file.write(header)
    lenNorm_file.write(header)
#   if psi_output_file:
#       psi_output_file.write(header)
    for event_key in event2sample2raw_counts:
        raw_counts_list = [event2sample2raw_counts[event_key][first_sample]]
        lenNorm_counts_list = [event2sample2lenNorm_counts[event_key][first_sample]]
        
#        psi_list = [getPSI(event2sample2counts[event_key][first_sample])]
        
        for samp in samples:
            if samp in event2sample2raw_counts[event_key]:
                raw_counts_list.append(event2sample2raw_counts[event_key][samp])
                lenNorm_counts_list.append(event2sample2lenNorm_counts[event_key][samp])
#                psi_list.append(getPSI(event2sample2counts[event_key][samp]))
            else:
                raw_counts_list.append("0;0")
                lenNorm_counts_list.append("0;0")
#                psi_list.append(NA)

        if "intron_retention" in event_key:
            if not hasInclusionCounts(event_key, 
                                      left_intron2sample2raw_counts,
                                      right_intron2sample2raw_counts):
                skipped_ir_events.add(event_key)
                continue

        out_key = getOutKey(event_key, 
                            event2genesStrand[event_key][0],
                            event2genesStrand[event_key][1],
                            event2genesStrand[event_key][2])

        if len(raw_counts_list) != num_samples:
            print "Error: Issue with raw counts for %s" % out_key
            sys.exit(1)
        if len(lenNorm_counts_list) != num_samples:
            print "Error: Issue with lenNorm counts for %s" % out_key
            sys.exit(1)
            

        raw_outline = "%s\t%s\n" % (out_key, 
                                "\t".join(raw_counts_list))
        lenNorm_outline = "%s\t%s\n" % (out_key, 
                                "\t".join(lenNorm_counts_list))

        output_file.write(raw_outline)
        lenNorm_file.write(lenNorm_outline)
#       if psi_output_file:
#           outline = "%s\t%s\n" % (out_key,
#                                   "\t".join(psi_list))
#           psi_output_file.write(outline)

    output_file.close() 
    lenNorm_file.close() 

#   if psi_output_file:
#       psi_output_file.close()

    # Printing left counts
    left_file = open(left_file_name, "w")
    lenNorm_left_file = open(lenNorm_left_file_name, "w")
    left_file.write(header)
    lenNorm_left_file.write(header)

#   psi_left_file = None
#   if psi_left_file_name:
#       psi_left_file = open(psi_left_file_name, "w")
#       psi_left_file.write(header)

    for event_key in left_intron2sample2raw_counts:
        if event_key in skipped_ir_events:
            continue

        raw_counts_list = [left_intron2sample2raw_counts[event_key][first_sample]]
        lenNorm_counts_list = [left_intron2sample2lenNorm_counts[event_key][first_sample]]
#        psi_list = [getPSI(left_intron2sample2counts[event_key][first_sample])]
        for samp in samples:
            if samp in left_intron2sample2raw_counts[event_key]:
                raw_counts_list.append(left_intron2sample2raw_counts[event_key][samp])
                lenNorm_counts_list.append(left_intron2sample2lenNorm_counts[event_key][samp])
#                psi_list.append(getPSI(left_intron2sample2counts[event_key][samp]))
            else:
                raw_counts_list.append("0;0")
                lenNorm_counts_list.append("0;0")
#                psi_list.append(NA)

        out_key = getOutKey(event_key, 
                            event2genesStrand[event_key][0],
                            event2genesStrand[event_key][1],
                            event2genesStrand[event_key][2])


        if len(raw_counts_list) != num_samples:
            print "Error: Issue with raw counts for %s" % out_key
            sys.exit(1)
        if len(lenNorm_counts_list) != num_samples:
            print "Error: Issue with lenNorm counts for %s" % out_key
            sys.exit(1)

        raw_outline = "%s\t%s\n" % (out_key,
                                "\t".join(raw_counts_list))
        lenNorm_outline = "%s\t%s\n" % (out_key,
                                "\t".join(lenNorm_counts_list))

        left_file.write(raw_outline)
        lenNorm_left_file.write(lenNorm_outline)

#       if psi_left_file:
#           outline = "%s\t%s\n" % (out_key,
#                                   "\t".join(psi_list))
#           psi_left_file.write(outline)

    left_file.close()
    lenNorm_left_file.close()

#   if psi_left_file:
#       psi_left_file.close()


    # Printing right counts
    right_file = open(right_file_name, "w")
    lenNorm_right_file = open(lenNorm_right_file_name, "w")
    right_file.write(header)
    lenNorm_right_file.write(header)

#   psi_right_file = None
#   if psi_right_file_name:
#       psi_right_file = open(psi_right_file_name, "w")
#       psi_right_file.write(header)

    for event_key in right_intron2sample2raw_counts:
        if event_key in skipped_ir_events:
            continue
        raw_counts_list = [right_intron2sample2raw_counts[event_key][first_sample]]
        lenNorm_counts_list = [right_intron2sample2lenNorm_counts[event_key][first_sample]]
#        psi_list = [getPSI(right_intron2sample2counts[event_key][first_sample])]
        for samp in samples:
            if samp in right_intron2sample2raw_counts[event_key]:
                raw_counts_list.append(right_intron2sample2raw_counts[event_key][samp])
                lenNorm_counts_list.append(right_intron2sample2lenNorm_counts[event_key][samp])
#                psi_list.append(getPSI(right_intron2sample2counts[event_key][samp]))
            else:
                raw_counts_list.append("0;0")
                raw_counts_list.append("0;0")
#                psi_list.append(NA)

        out_key = getOutKey(event_key, 
                            event2genesStrand[event_key][0],
                            event2genesStrand[event_key][1],
                            event2genesStrand[event_key][2])

        if len(raw_counts_list) != num_samples:
            print "Error: Issue with raw counts for %s" % out_key
            sys.exit(1)
        if len(lenNorm_counts_list) != num_samples:
            print "Error: Issue with lenNorm counts for %s" % out_key

        raw_outline = "%s\t%s\n" % (out_key,
                                "\t".join(raw_counts_list))
        lenNorm_outline = "%s\t%s\n" % (out_key,
                                "\t".join(lenNorm_counts_list))

        right_file.write(raw_outline)
        lenNorm_right_file.write(lenNorm_outline)

#       if psi_right_file:
#           outline = "%s\t%s\n" % (out_key,
#                                   "\t".join(psi_list))
#           psi_right_file.write(outline)

    right_file.close()
    lenNorm_right_file.close()

#   if psi_right_file:
#       psi_right_file.close()
     
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

def getKeyandCount(line):
    line_list = line.split("\t")

#    event_key_list = line_list[0:1] + line_list[2:3] + line_list[4:5] + line_list[6:11]
    event_key_list = line_list[0:2] + line_list[3:4] + line_list[5:10]
    event_key = "\t".join(event_key_list)

    genes = line_list[2].split(",")
    strand = line_list[4]
    const_region = line_list[10]

    excl_raw_count = line_list[-4]
    incl_raw_count = line_list[-3] 
    excl_lenNorm_count = line_list[-2]
    incl_lenNorm_count = line_list[-1] 
    
    raw_count_str = "%s;%s" % (excl_raw_count, incl_raw_count)
    lenNorm_count_str = "%s;%s" % (excl_lenNorm_count, incl_lenNorm_count)

    if (("-" in raw_count_str) or ("-" in lenNorm_count_str)):
        print "Negative value in %s" % line
        raw_count_str = "0;0"
        lenNorm_count_str = "0;0"

    return event_key, raw_count_str, lenNorm_count_str, genes, strand, const_region

def getIRKeyandCounts(line, jcn_seq_len, lengthNorm):
    """ 
    Returns excl;incl pairs for the left side and then the right side
    """

    line_list = line.split("\t")

    if lengthNorm:
        excl_count = line_list[-2]
    else:
        excl_count = line_list[-4]

    left_right_incl_cts = line_list[19] 

    left_incl_count, right_incl_count = left_right_incl_cts.split(";")

    if lengthNorm:
        left_incl_count = normalizeByLen(int(left_incl_count), jcn_seq_len)
        right_incl_count = normalizeByLen(int(right_incl_count), jcn_seq_len)
    
        left_count_str = "%s;%d" % (excl_count, left_incl_count)
        right_count_str = "%s;%d" % (excl_count, right_incl_count)
    else:
        left_count_str = "%s;%s" % (excl_count, left_incl_count)
        right_count_str = "%s;%s" % (excl_count, right_incl_count)

    return left_count_str, right_count_str

def getOutKey(event_key, genes, strand, const_region):
    """
    Inserts the gene and strand information into the output events
    """
    event_list = event_key.split("\t") 
    gene_names = list(genes)
    gene_names.sort()
    event_list.insert(2, ",".join(gene_names))
    event_list.insert(4, strand)
    event_list.append(const_region)

    return "\t".join(event_list)

def getPSI(excl_incl_ct_str):

    excl_str, incl_str = excl_incl_ct_str.split(";")

    try:
        excl = float(excl_str)
        incl = float(incl_str)
    except:
        print "Warning:Bad PSI value"
        return NA

    if excl + incl == 0:
        return NA

    psi = (incl/(incl + excl)) * 100

    psi_str = "%.2f" % psi

    return psi_str

def hasInclusionCounts(event_key,
                       left_intron2sample2counts,
                       right_intron2sample2counts):
    leftHasInclusion = False
    rightHasInclusion = False

    for samp in left_intron2sample2counts[event_key]:
        count_str = left_intron2sample2counts[event_key][samp]
        excl_ct, incl_ct = map(int,count_str.split(";"))

        if incl_ct > 0:
            leftHasInclusion = True
            break

    for samp in right_intron2sample2counts[event_key]:
        count_str = right_intron2sample2counts[event_key][samp]
        excl_ct, incl_ct = map(int,count_str.split(";"))

        if incl_ct > 0:
            rightHasInclusion = True
            break

    return (leftHasInclusion and rightHasInclusion)

def updateGeneStrand(event2genesStrand, event_key, genes, strand, const_region):
    """
    Will resolve strand and gene information.
    Dictionary is of the format:
    {event_key: (set([genes,]),
                 strand,
                 const_region)}
    """
    try:
        this_strand = event2genesStrand[event_key][1]
        if this_strand != strand:
            event2genesStrand[event_key][1] = "."

        this_const_region = event2genesStrand[event_key][2]
        if this_const_region != const_region:
            if this_const_region == "":
                event2genesStrand[event_key][2] = const_region
    except:
        # Key does not exist
        event2genesStrand[event_key] = (set([]),
                                        strand,
                                        const_region)

    genes2update = []
    for gene in genes:
        if gene != "None":
            genes2update.append(gene)

    event2genesStrand[event_key][0].update(genes2update)

    
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
