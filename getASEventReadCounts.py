#!/lab/64/bin/python
# getASEventReadCounts.py
# Author: Angela Brooks
# Program Completion Date:
# Description: Takes junction read counts between two samples and outputs:
# A file of cassette events:
# K/N (known vs. novel) chr exon_start exon_end exclusion_count1
# inclusion_count1 exclusion_count2 inclusion_count2
# A file of alternative donor events:
# K/N (known vs. novel) chr included_donor_site_intron
# excluded_donor_site_intron(s) exclusion_count1 inclusion_count1
# exclusion_count2 inclusion_count2
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

from __future__ import division

import sys
import optparse 
import pdb
import os
import profile
import pickle

from helperFunctions import updateDictOfLists, updateDictOfSets, coordsOverlap, runCmd
from coord_helperFunctions import getSearchTree, hasOuterContainer, findInternalCoords, hasOverlap

from Bio import SeqIO
#############
# CONSTANTS #
#############
HOST = "localhost"
USER = "root"
PASSWD = ""

# This should be a float for easier use in division
DEF_EXON_LEN_NORM = 100.0
DO_LEN_NORM = True

INFINITY = 10000000000000000000000000000000000000
NOVEL_DIST = 2000

# Based on longest exon size in 
MAX_EXON_LEN = 30000
MAX_GENE_LEN = 400000

# For mutually exclusive and multi-cassette exons, the calculation explodes
# when there are too many possible exons in the cluster.  For now, I am
# filtering out any events that have a large number of exons.
MAX_EXON_CLUSTER = 12

READ_LENGTH = 37

# FISHER'S TEST SCRIP
BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
FISHER_SCRIPT = "%s/fisherTest_exclusionInclusion.py" % BIN_DIR
if not os.path.exists(FISHER_SCRIPT):
    print "ERROR: fisherTest_exclusionInclusion.py needs to be in the same directory."
    sys.exit(1)
IR_SCRIPT = "%s/getIntronRetentionEvents.py" % BIN_DIR
if not os.path.exists(IR_SCRIPT):
    print "ERROR: getIntronRetentionEvents.py needs to be in the same directory."
    sys.exit(1)
COORD_COUNT = "%s/coordReadCounts.py" % BIN_DIR
if not os.path.exists(COORD_COUNT):
    print "ERROR: coordReadCounts.py needs to be in the same directory."
    sys.exit(1)

SHELL = "/bin/tcsh"

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
    opt_parser.add_option("--jcn1",
                          dest="jcn1",
                          type="string",
                          help="""First file of jcn counts (e.g. untreated)
                                  From single read data. BED format.""",
                          default=None)
    opt_parser.add_option("--jcn2",
                          dest="jcn2",
                          type="string",
                          help="""Second file of jcn counts (e.g. RNAi).  From
                                  single read data. BED format""",
                          default=None)
    opt_parser.add_option("--genome_reads1",
                          dest="genome_reads1",
                          type="string",
                          help="""Alignment output of reads to genome only from
                                  the first sample""",
                          default=None)
    opt_parser.add_option("--genome_reads2",
                          dest="genome_reads2",
                          type="string",
                          help="""Alignment output of reads to genome only from
                                  the first sample""",
                          default=None)
    opt_parser.add_option("--coord_counts1",
                          dest="coord_counts1",
                          type="string",
                          help="""Instead of giving reads as input, can give
                                  precomputed output from coordReadCounts.py    
                                  from the first sample""",
                          default=None)
    opt_parser.add_option("--coord_counts2",
                          dest="coord_counts2",
                          type="string",
                          help="""Instead of giving reads as input, can give
                                  precomputed output from coordReadCounts.py    
                                  from the first sample""",
                          default=None)
    opt_parser.add_option("--ie1",
                          dest="intron_exon1",
                          type="string",
                          help="""First file (e.g. untreated) of counts to exon/intron
                                junctions.  The first column will be the intron
                                coordinate.  The second column are the counts to
                                the left side of the junction.  The third
                                column are the counts to the right side of the
                                junction""",
                          default=None)
    opt_parser.add_option("--ie2",
                          dest="intron_exon2",
                          type="string",
                          help="Option: Second file of exon/intron counts (e.g. RNAi)",
                          default=None)
    opt_parser.add_option("--sqlite_db_dir",
                          dest="sqlite_db_dir",
                          type="string",
                          help="""Location of sqlite databases. If sqlite
                                  databases are used, will override usage of a 
                                  MySQL database.""",
                          default=None)
    opt_parser.add_option("--host",
                          dest="host",
                          type="string",
                          help="MySQL database host. Def=\'%s\'" % HOST,
                          default=HOST)
    opt_parser.add_option("--user",
                          dest="user",
                          type="string",
                          help="MySQL database user. Def=\'%s\'" % USER,
                          default=USER)
    opt_parser.add_option("--passwd",
                          dest="passwd",
                          type="string",
                          help="MySQL database password. Def=\'%s\'" % PASSWD,
                          default=PASSWD)
    opt_parser.add_option("--txt_db1",
                          dest="txt_db1",
                          type="string",
                          help="""Database of transcript annotations derived
                                  from a gtf file. Used to define exon and
                                  intron and gene coordinates.""",
                          default=None)
    opt_parser.add_option("--txt_db2",
                          dest="txt_db2",
                          type="string",
                          help="""Database of transcript annotations derived
                                  from a gtf file. Used to identify alternative
                                  first and last exons.  Can be the same or
                                  different as txt_db1.  This annotation should
                                  be fairly clean of fragmented
                                  transcripts.""",
                          default=None)
    opt_parser.add_option("--txt_db3",
                          dest="txt_db3",
                          type="string",
                          help="""Database of transcript annotations derived
                                  from a gtf file. Used for annotating gene
                                  names and whether an intron/junction is annotated or
                                  not. By default, txt_db1 will be used for this
                                  information.""",
                          default=None)
    opt_parser.add_option("-p",
                          dest="prefix",
                          type="string",
                          help="Prefix string to output files.",
                          default=None)
    opt_parser.add_option("--norm1",
                          dest="norm1",
                          type="float",
                          help="""Normalization factor to divide the first sample
                                  counts by.  Typically, this would be the
                                  total number of reads in the sample.""",
                          default=None)
    opt_parser.add_option("--norm2",
                          dest="norm2",
                          type="float",
                          help="""Normalization factor to divide the second  sample
                                  counts by.  Typically, this would be the
                                  total number of reads in the sample.""",
                          default=None)
#   opt_parser.add_option("--lengthNorm",
#                         dest="lengthNorm",
#                         action="store_true",
#                         help="""Default is to not normalize read counts by
#                                isoform length. This will option will specify
#                                normalize read counts by isoform length.""",
#                         default=False)
    opt_parser.add_option("--jcn_seq_len",
                          dest="jcn_seq_len",
                          type="int",
                          help="""I recommmend this value to be
                                  (read_length-6)*2 which assumes reads aligned
                                  to junctions with at least a 6pb overhang. If
                                  tophat was used for the alignment and you
                                  used the -a option with something < 6, give
                                  the value (read_length-(anchor length)*2.""",
                          default=None)
#   opt_parser.add_option("--method",
#                         dest="method",
#                         type="string",
#                         help="""Type of correction method:
#                                 'BH' - Benjamini & Hochberg,
#                                 'bonferroni'""",
#                         default=None)
    opt_parser.add_option("--fasta",
                           dest="genome_file",
                           type="string",
                           help="""Contains the genome sequence organized by
                                   chromosome.""",
                           default=None)
    opt_parser.add_option("--disambiguate_jcn_strand",
                           dest="disambiguate_jcn_strand",
                           action="store_true",
                           help="""For junction alignments that had an
                                   undetermined strand, this flag will
                                   determine the strand of the junction based on
                                   the genome sequence at the splice sites.""",
                           default=False)
    opt_parser.add_option("--by_chr",
                           dest="this_chr",
                           type="string",
                           help="""If running the script one chromosome at a
                                   time, this will specify which chromosome.""",
                           default=False)
    opt_parser.add_option("--keep_intermediate",
                           dest="keep_interm",
                           action="store_true",
                           help="""Will remove intermediate files by default. Use
                                    this option to keep them.""",
                           default=False)

    # Checking paths of constants
    if not os.path.exists(FISHER_SCRIPT):
        print "Change location of the fisher's script: %s" % FISHER_SCRIPT
        opt_parser.print_help()
        sys.exit(1)
    if not os.path.exists(IR_SCRIPT):
        print "Change location of the IR script: %s" % IR_SCRIPT
        opt_parser.print_help()
        sys.exit(1)
    if not os.path.exists(COORD_COUNT):
        print "Change location of coordReadCount.py: %s" % COORD_COUNT
        opt_parser.print_help()
        sys.exit(1)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--jcn1")
    opt_parser.check_required("--jcn2")
    opt_parser.check_required("--txt_db1")
    opt_parser.check_required("--txt_db2")
#    opt_parser.check_required("--method")
    opt_parser.check_required("--jcn_seq_len")

#    method = options.method

    jcn_seq_len = options.jcn_seq_len
#   if options.lengthNorm:
#       global DO_LEN_NORM
#       DO_LEN_NORM = True

#   if method != "BH" and method != "bonferroni":
#       print "Wrong method given."
#       opt_parser.print_help()
#       sys.exit(1)

    keep_interm = options.keep_interm

    genome_file = options.genome_file
    disambiguate_jcn_strand = options.disambiguate_jcn_strand

    this_chr = None
    if options.this_chr:
        this_chr = options.this_chr

    jcn1_file = None
    jcn2_file = None
    if options.jcn1:
        jcn1_file = open(options.jcn1)
        jcn2_file = open(options.jcn2)

    countExonReads = False
    printExonCoords = False
    
    genome_read_file1 = None
    genome_read_file2 = None
    coord_counts1 = None
    coord_counts2 = None

    exon_coords = None
    if options.genome_reads1 or options.coord_counts1:
        countExonReads = True
        printExonCoords = True
        exon_coords = set([])

    if options.genome_reads1:
        if options.coord_counts1:
            print "Must select either genome reads as input or coord counts."
            opt_parser.print_help()
            sys.exit(1)
        genome_read_file1 = options.genome_reads1
    elif options.coord_counts1:
        coord_counts1 = options.coord_counts1

    if options.genome_reads2:
        if options.coord_counts2:
            print "Must select either genome reads as input or coord counts."
            opt_parser.print_help()
            sys.exit(1)
        genome_read_file2 = options.genome_reads2
    elif options.coord_counts2:
        coord_counts2 = options.coord_counts2

    paired_genome_read_file1 = None
    paired_genome_read_file2 = None
    paired_coord_counts1 = None
    paired_coord_counts2 = None
    paired_read_w_coord1 = None
    paired_read_w_coord2 = None

    paired_junctions2qname_file1 = None
    paired_junctions2qname_file2 = None
    paired_ie_junctions2qname_file1 = None
    paired_ie_junctions2qname_file2 = None
    
    paired_coord2qname2count1 = None
    paired_coord2qname2count2 = None
    paired_junction2qname2count1 = None
    paired_junction2qname2count2 = None
    paired_ie_junction2qname2count1 = None
    paired_ie_junction2qname2count2 = None

    paired_endCounting1 = False
    paired_endCounting2 = False

    ie1_file = None
    ie2_file = None
    if options.intron_exon2 is not None and options.intron_exon2 is not None:
        ie1_file = open(options.intron_exon1)
        ie2_file = open(options.intron_exon2)

    # Remnant of old command line options
    doFisher = False
    a_exons_only = True

    txt_db1 = options.txt_db1
    txt_db2 = options.txt_db2

    annot_db = txt_db1
    if options.txt_db3:
        annot_db = options.txt_db3

    db = None
    if options.sqlite_db_dir:
        import pysqlite_wrap
        db = pysqlite_wrap.DB(formatDir(options.sqlite_db_dir))
    else: # Use MySQL database
        import mysqldb_wrap
        db = mysqldb_wrap.DB(options.host, options.user, options.passwd)

#    read_length = options.read_length

    norm1 = options.norm1
    norm2 = options.norm2

    prefix = None
    if options.prefix:
        prefix = options.prefix
        # Output files
        cassette_out = open(prefix + "_cassette_event_counts.txt", "w")
        donor_out = open(prefix + "_alternative_donor_counts.txt", "w")
        accept_out = open(prefix + "_alternative_acceptor_counts.txt", "w")
        jcn_only_donor_out = open(prefix + "_jcn_only_alternative_donor_counts.txt", "w")
        jcn_only_accept_out = open(prefix + "_jcn_only_alternative_acceptor_counts.txt", "w")
        afe_out = open(prefix + "_alternative_first_exon_counts.txt", "w")
        ale_out = open(prefix + "_alternative_last_exon_counts.txt", "w")
        me_out = open(prefix + "_mutually_exclusive_counts.txt", "w")
        mc_out = open(prefix + "_coord_cassette_counts.txt", "w")
#       if countExonReads:
#           polya_out = open(prefix + "_alternative_polyA_counts.txt", "w")
        if ie1_file is not None:
            ir_left_out = open(prefix + "_intron_retention_left_counts.txt", "w")
            ir_right_out = open(prefix + "_intron_retention_right_counts.txt", "w")

        all_event_info_out = open(prefix + "_all_AS_event_info.txt", "w")

        if __name__ == "__main__":
            ERROR_LOG = open("%s_error.log" % prefix, "w")
    else:
        # Output files
        cassette_out = open("cassette_event_counts.txt", "w")
        donor_out = open("alternative_donor_counts.txt", "w")
        accept_out = open("alternative_acceptor_counts.txt", "w")
        jcn_only_donor_out = open("jcn_only_alternative_donor_counts.txt", "w")
        jcn_only_accept_out = open("jcn_only_alternative_acceptor_counts.txt", "w")
        afe_out = open("alternative_first_exon_counts.txt", "w")        
        ale_out = open("alternative_last_exon_counts.txt", "w")        
        me_out = open("mutually_exclusive_counts.txt", "w")        
        mc_out = open("coord_cassette_counts.txt", "w")        
#       if countExonReads:
#           polya_out = open("alternative_polyA_counts.txt", "w")        
        if ie1_file is not None:
            ir_left_out = open("intron_retention_left_counts.txt", "w")
            ir_right_out = open("intron_retention_right_counts.txt", "w")

        all_event_info_out = open("all_AS_event_info.txt", "w")

        if __name__ == "__main__":
            ERROR_LOG = open("error.log", "w")
    # jcn_count_dict - {jcn_coord_str:(file1_count, file2_count)}
    # coord_start2end - {chr:{start:[end]}}
    # coord_end2start - {chr:{end:[start]}}
    print "Parsing input."

    jcn_count_dict = None
    coord_start2end = None
    coord_end2start = None
    
    # Search tree: {chr:strand:searchTreeOfCoords}
    (all_jcn_count_dict,
     all_coord_start2end,
     all_coord_end2start,
     all_jcn2strand,
     all_jcn_search_tree) = parseJcns(jcn1_file, jcn2_file, genome_file, disambiguate_jcn_strand)

    jcn_chrs = all_jcn_search_tree.keys()

    # full_exon_count_dict = {chr_start_end:count}
    # start_exon_count_dict = {chr_start:sum of counts}
    # end_exon_count_dict = {chr_end:sum of counts}
    full_exon_count_dict = None
    start_exon_count_dict = None
    end_exon_count_dict = None   
    full_multi_exon_count_dict = None
    start_multi_exon_count_dict = None
    end_multi_exon_count_dict = None

    ir_count_dict = None
    if ie1_file is not None:
        # ir_count_dict - {jcn_coord_str: {left: (file1_count, file2_count),
        #                                  right: {file1_count, file2_count)}
        ir_count_dict = parseIEJcnFiles(ie1_file, ie2_file)

    # {chr: set[(start, end, strand)])
    annotated_introns = getAnnotatedIntronCoords(db, annot_db, this_chr)

    chrCheck(jcn_chrs, annotated_introns, "txt_db1:annotated introns")

    annotated_exons = None
    if a_exons_only:
        # {chr: set[(start, end, strand)])
        (annotated_exons,
         annotated_internal_exons,
         annotated_exons_no_strand,
         annotated_exons_by_strand,
         annotated_exon_search_tree) = getAnnotatedExonCoords(db, txt_db1, this_chr)

    chrCheck(jcn_chrs, annotated_exons, "txt_db1:annotated exons")
    chrCheck(jcn_chrs, annotated_internal_exons, None)
    chrCheck(jcn_chrs, annotated_exons_no_strand, None)
    chrCheck(jcn_chrs, annotated_exons_by_strand, None)
    chrCheck(jcn_chrs, annotated_exon_search_tree, None)

    # {chr: (start,end): [gene_names]}
    # {chr: strand: (start,end): [gene_names]}
    (annotated_genes,
     annotated_genes_by_strand) = getAnnotatedGenes(db, annot_db, this_chr)

    chrCheck(jcn_chrs, annotated_genes, "txt_db1:annotated genes")
    chrCheck(jcn_chrs, annotated_genes_by_strand, None)

    # {chr:(start, end)} 
    (alt_first_exons, 
     alt_last_exons, alt_first_exons_start2end,
     alt_first_exons_end2start,
     alt_last_exons_start2end,
     alt_last_exons_end2start,
     alt_first_exon_search_tree,
     alt_last_exon_search_tree) = getFirstLastExons(db, txt_db2, this_chr)

    if not annotated_exons:
        print "ERROR: No exon entries in --txt_db1. Please check the \"exon\" table"
        sys.exit(1)

    if annotated_genes == {}:
        print "ERROR: No genes in --txt_db1. Please check the \"gene\" table"
        sys.exit(1)

    if alt_first_exons == {} or alt_last_exons == {}:
        print "ERROR: Error with transcripts in --txt_db2. Please check this database."
        sys.exit(1)

    chrCheck(jcn_chrs, alt_first_exons, "txt_db2:alt first exons")
    chrCheck(jcn_chrs, alt_last_exons, "txt_db2:alt last exons")
    chrCheck(jcn_chrs, alt_first_exons_start2end, None)
    chrCheck(jcn_chrs, alt_first_exons_end2start, None)
    chrCheck(jcn_chrs, alt_last_exons_start2end, None)
    chrCheck(jcn_chrs, alt_last_exons_end2start, None)
    chrCheck(jcn_chrs, alt_first_exon_search_tree, None)
    chrCheck(jcn_chrs, alt_last_exon_search_tree, None)
        

    # A merge of internal exons from txt_db1 and alternative first and last
    # exons from txt_db2
    all_confident_exons = {}
    all_confident_exons_start2end = {}
    all_confident_exons_end2start = {}
    for chr in annotated_internal_exons:
        all_confident_exons[chr] = set([])
        all_confident_exons[chr].update(set(annotated_internal_exons[chr]))

        all_confident_exons_start2end[chr] = {}
        all_confident_exons_end2start[chr] = {}

        for (start, end) in annotated_internal_exons[chr]:
            updateDictOfSets(all_confident_exons_start2end[chr], start, end)
            updateDictOfSets(all_confident_exons_end2start[chr], end, start)

        try: # In case the chr do not exist in alt_first or last exons
            all_confident_exons[chr].update(set(alt_first_exons[chr]))
            
            for (start, end) in alt_first_exons[chr]:
                updateDictOfSets(all_confident_exons_start2end[chr], start, end)
                updateDictOfSets(all_confident_exons_end2start[chr], end, start)
        except:
            continue
        try:
            all_confident_exons[chr].update(set(alt_last_exons[chr]))

            for (start, end) in alt_last_exons[chr]:
                updateDictOfSets(all_confident_exons_start2end[chr], start, end)
                updateDictOfSets(all_confident_exons_end2start[chr], end, start)
        except:
            continue

    # In case alt first exons not on a chr in internal exons
    for chr in alt_first_exons:
        if chr not in all_confident_exons:
            all_confident_exons[chr] = set([]) 
            all_confident_exons[chr].update(set(alt_first_exons[chr]))

            all_confident_exons_start2end[chr] = {}
            all_confident_exons_end2start[chr] = {}

            for (start, end) in alt_first_exons[chr]:
                updateDictOfSets(all_confident_exons_start2end[chr], start, end)
                updateDictOfSets(all_confident_exons_end2start[chr], end, start)
            try:
                all_confident_exons[chr].update(set(alt_last_exons[chr]))
            except:
                continue
    # In case alt last exons not on a chr in internal exons
    for chr in alt_last_exons:
        if chr not in all_confident_exons:
            all_confident_exons[chr] = set([]) 
            all_confident_exons[chr].update(set(alt_last_exons[chr]))

            all_confident_exons_start2end[chr] = {}
            all_confident_exons_end2start[chr] = {}

            for (start, end) in alt_last_exons[chr]:
                updateDictOfSets(all_confident_exons_start2end[chr], start, end)
                updateDictOfSets(all_confident_exons_end2start[chr], end, start)

    # Used to filter junction out of the altDonor/Acceptor events
    # These are all the junctions that skip exons.
    # excl_jcns = {chr:set(start, end)}
    excl_jcns = {}

    print "Cassette Exons"
    # cassette_exons = {chr:set(start,end)}
    cassette_exons = printCassetteExons(db, 
                                        alt_first_exons,
                                        alt_last_exons,
                                        annotated_genes,
                                        annotated_genes_by_strand,
                                        all_jcn_count_dict, 
                                        all_coord_start2end,
                                        all_coord_end2start,
                                        all_jcn2strand,
                                        full_exon_count_dict,
                                        full_multi_exon_count_dict,
                                        start_multi_exon_count_dict,
                                        end_multi_exon_count_dict,
                                        annotated_exons,
                                        annotated_exons_by_strand,
                                        annotated_introns,
                                        printExonCoords,
                                        exon_coords,
                                        excl_jcns,
                                        cassette_out,
                                        all_event_info_out,
                                        norm1, norm2,
                                        jcn_seq_len)

    cassette_out.close()

    print "Mutually Exclusive"
    # Cassette exons is updated in this function
    printMutuallyExclusive(db,
                           annotated_genes,
                           annotated_genes_by_strand,
                           annotated_internal_exons,
                           alt_first_exons,
                           alt_last_exons,
                           all_jcn_count_dict,
                           all_coord_start2end,
                           all_coord_end2start,
                           all_jcn2strand,
                           full_exon_count_dict,
                           full_multi_exon_count_dict,
                           cassette_exons,
                           annotated_introns,
                           annotated_exons,
                           annotated_exons_by_strand,
                           annotated_exon_search_tree,
                           printExonCoords,
                           exon_coords,
                           excl_jcns,
                           me_out,
                           all_event_info_out,
                           norm1, norm2, jcn_seq_len)

    me_out.close()

    print "Multi-Cassette"
    printMultiCassetteExons(db,
                            annotated_genes,
                            annotated_genes_by_strand,
                            alt_first_exons,
                            alt_last_exons,
                            alt_first_exon_search_tree,
                            alt_last_exon_search_tree,
                            all_jcn_count_dict, 
                            all_coord_start2end,
                            all_coord_end2start,
                            all_jcn2strand,
                            all_jcn_search_tree,
                            full_exon_count_dict,
                            full_multi_exon_count_dict,
                            cassette_exons,
                            annotated_introns,
                            annotated_exons,
                            annotated_exons_by_strand,
                            printExonCoords,
                            exon_coords,
                            excl_jcns,
                            mc_out,
                            all_event_info_out,
                            norm1, norm2, jcn_seq_len)

    mc_out.close()

    print "Alternative Donors and Acceptors"
    printAlternativeDonorsAcceptors(db,
                                    annotated_genes,
                                    annotated_genes_by_strand,
                                    alt_first_exons,
                                    alt_last_exons,
                                    alt_first_exons_start2end,
                                    alt_first_exons_end2start,
                                    alt_last_exons_start2end,
                                    alt_last_exons_end2start,
                                    all_jcn_count_dict,
                                    all_coord_start2end,
                                    all_coord_end2start,
                                    all_jcn2strand,
                                    ir_count_dict,
                                    cassette_exons,
                                    annotated_introns,
                                    annotated_exons_no_strand,
                                    annotated_exons_by_strand,
                                    all_confident_exons,
                                    all_confident_exons_start2end,
                                    all_confident_exons_end2start,
                                    printExonCoords,
                                    exon_coords,
                                    excl_jcns,
                                    donor_out, afe_out, jcn_only_donor_out, 
                                    accept_out, ale_out, jcn_only_accept_out,
                                    all_event_info_out,
                                    norm1, norm2)




    donor_out.close()
    afe_out.close()
    jcn_only_donor_out.close()
    accept_out.close()
    ale_out.close()
    jcn_only_accept_out.close()

    print "Intron Retention Events"
    if ir_count_dict is not None:
        printIREvents(db, annotated_genes, annotated_genes_by_strand, annotated_exons,
                      annotated_exons_by_strand, all_coord_start2end, all_coord_end2start,
                      all_jcn_count_dict, all_jcn2strand, ir_count_dict, 
                      ir_left_out, ir_right_out, printExonCoords, exon_coords,
                      norm1, norm2, jcn_seq_len)

    if ie1_file is not None:
        ir_left_out.close()
        ir_right_out.close()

    all_event_info_out.close()

    if options.prefix:
        prefix = options.prefix
        # Output files
        cassette_out_str = prefix + "_cassette_event_counts.txt"
        donor_out_str = prefix + "_alternative_donor_counts.txt"
        accept_out_str = prefix + "_alternative_acceptor_counts.txt"
        jcn_only_donor_out_str = prefix + "_jcn_only_alternative_donor_counts.txt"
        jcn_only_accept_out_str = prefix + "_jcn_only_alternative_acceptor_counts.txt"
        afe_out_str = prefix + "_alternative_first_exon_counts.txt"
        ale_out_str = prefix + "_alternative_last_exon_counts.txt"
        me_out_str = prefix + "_mutually_exclusive_counts.txt"
        mc_out_str = prefix + "_coord_cassette_counts.txt"

        if ie1_file is not None:
            ir_left_out_str = prefix + "_intron_retention_left_counts.txt"
            ir_right_out_str = prefix + "_intron_retention_right_counts.txt"

        all_event_info_str = prefix + "_all_AS_event_info.txt"
    else:
        # Output files
        cassette_out_str = "cassette_event_counts.txt"
        donor_out_str = "alternative_donor_counts.txt"
        accept_out_str = "alternative_acceptor_counts.txt"
        jcn_only_donor_out_str = "jcn_only_alternative_donor_counts.txt"
        jcn_only_accept_out_str = "jcn_only_alternative_acceptor_counts.txt"
        afe_out_str = "alternative_first_exon_counts.txt"
        ale_out_str = "alternative_last_exon_counts.txt"
        me_out_str = "mutually_exclusive_counts.txt"
        mc_out_str = "coord_cassette_counts.txt"

        if ie1_file is not None:
            ir_left_out_str = "intron_retention_left_counts.txt"
            ir_right_out_str = "intron_retention_right_counts.txt"

        all_event_info_str = "all_AS_event_info.txt"

    if printExonCoords:
        if genome_read_file1 or genome_read_file2:
            # Produce exon coord file
            if prefix:            
                exon_coord_file_name = "%s_tmp_exon_coord_file.txt" % prefix
            else:
                exon_coord_file_name = "tmp_exon_coord_file.txt"
    
            exon_coord_file = open(exon_coord_file_name, "w")

            for coord in exon_coords:

                s = "%s\t%d\t%d\n" % (coord[0],
                                    coord[1],
                                    coord[2])
                exon_coord_file.write(s)

            exon_coord_file.close()

        if genome_read_file1:
            if prefix:
                mapped_file1_name = "%s_tmp_exon_coord_file1_readCounts.txt"  % prefix
            else:
                mapped_file1_name = "tmp_exon_coord_file1_readCounts.txt" 

            cmd = "python %s --coords %s " % (COORD_COUNT, exon_coord_file_name)
            
            # Run First Set of Reads through coordReadCounts
            first_cmd = cmd + "--reads %s -o %s" % (genome_read_file1,
                                                    mapped_file1_name) 
                                
            print "Running: %s" % first_cmd
            os.system(first_cmd)
        else:
            mapped_file1_name = coord_counts1

        if genome_read_file2:
            if prefix:
                mapped_file2_name = "%s_tmp_exon_coord_file2_readCounts.txt"  % prefix
            else:
                mapped_file2_name = "tmp_exon_coord_file2_readCounts.txt" 

            cmd = "python %s --coords %s " % (COORD_COUNT, exon_coord_file_name)

            # Run Second Set of Reads through coordReadCounts
            second_cmd = cmd + "--reads %s -o %s" % (genome_read_file2,
                                                     mapped_file2_name) 
                                
            print "Running: %s" % second_cmd
            os.system(second_cmd)

        else:
            mapped_file2_name = coord_counts2

        # {coord_str: count}
        mapped_file1_counts = parseCoordCounts(mapped_file1_name, norm1)
        mapped_file2_counts = parseCoordCounts(mapped_file2_name, norm2)

        # Counts calculated here need to be maintained
        # {exon_str: (excl1, incl1, excl2, incl2)}
        ce2total_counts = updateCounts2Cassette(cassette_out_str, 
                              mapped_file1_counts, mapped_file2_counts,
                              norm1, norm2, jcn_seq_len)

        # Will return a dictionary:
        # {(incl_str, excl_str): (excl1, incl1, excl2, incl2)}
        alt_donor2total_counts = updateCounts2AltDonorAccept(donor_out_str, 
                                    ir_count_dict,
                                    mapped_file1_counts,mapped_file2_counts,
                                    norm1, norm2, jcn_seq_len)
                                    
        alt_accept2total_counts = updateCounts2AltDonorAccept(accept_out_str, 
                                    ir_count_dict,
                                    mapped_file1_counts,mapped_file2_counts,
                                    norm1, norm2, jcn_seq_len)

        updateCounts2AltDonorAccept(jcn_only_donor_out_str,
                                    ir_count_dict,
                                    mapped_file1_counts,mapped_file2_counts,
                                    norm1, norm2, jcn_seq_len, True)


        updateCounts2AltDonorAccept(jcn_only_accept_out_str,
                                    ir_count_dict,
                                    mapped_file1_counts,mapped_file2_counts,
                                    norm1, norm2, jcn_seq_len, True)

        # {(incl_str, excl_str): (excl1, incl1, excl2, incl2)}
        afe2total_counts = updateCounts2AFE_ALE(afe_out_str, 
                                                all_jcn_count_dict,
                                                all_coord_start2end,
                                                all_coord_end2start,
                             mapped_file1_counts, mapped_file2_counts,
                             norm1, norm2, jcn_seq_len)

        # {(incl_str, excl_str): (excl1, incl1, excl2, incl2)}
        ale2total_counts = updateCounts2AFE_ALE(ale_out_str, 
                                                all_jcn_count_dict,
                                                all_coord_start2end,
                                                all_coord_end2start,
                             mapped_file1_counts, mapped_file2_counts,
                             norm1, norm2, jcn_seq_len)

        # Will return a dictionary:
        # {(incl_exon, excl_exons): (excl1, incl1, excl2, incl2)}
        mxe2total_counts = updateCounts2MutuallyExclusive(me_out_str, 
                                       all_jcn_count_dict,
                                       all_coord_start2end,
                                       all_coord_end2start,
                                       mapped_file1_counts, mapped_file2_counts,
                                       norm1, norm2, jcn_seq_len)

        # Will return a dictionary:
        # {(incl_exons, excl_jcn): (excl1, incl1, excl2, incl2)}
        mc2total_counts = updateCounts2MultiCassette(mc_out_str, 
                                   all_jcn_count_dict,
                                   all_coord_start2end,
                                   all_coord_end2start,
                                   mapped_file1_counts, mapped_file2_counts,
                                   norm1, norm2, jcn_seq_len)

        updateCounts2all_as_events(all_event_info_str, 
                                   mapped_file1_counts, mapped_file2_counts,
                                   norm1, norm2, jcn_seq_len)

    # Now do intron retention events
    if ie1_file is not None:

        # Pickle up existing annoation dictionaries to make things more
        # efficient
        if prefix:
            annot_intron_pick_name = "%s_annotated_introns.pk" % prefix
        else:
            annot_intron_pick_name = "annotated_introns.pk"
        intron_pick_file = open(annot_intron_pick_name, "w")
        pickle.dump(annotated_introns, intron_pick_file)
        intron_pick_file.close()

        if prefix:
            annot_exon_pick_name = "%s_annoated_exons_by_strand.pk" % prefix
        else:
            annot_exon_pick_name = "annoated_exons_by_strand.pk"
        exon_pick_file = open(annot_exon_pick_name, "w")
        pickle.dump(annotated_exons_by_strand, exon_pick_file)
        exon_pick_file.close()

        if options.prefix:
            op = "-l %s_intron_retention_left_counts.txt " % prefix
            op += "-r %s_intron_retention_right_counts.txt " % prefix
            op += "-p %s " % prefix
            op += "--all_as_event %s_all_AS_event_info_irOnly.txt" % prefix
        else:
            op = "-l intron_retention_left_counts.txt "
            op += "-r intron_retention_right_counts.txt "
            op += "--all_as_event all_AS_event_info_irOnly.txt"

        cmd = "python %s %s " % (IR_SCRIPT, op)
        cmd += "--intron %s --exon %s" % (annot_intron_pick_name, annot_exon_pick_name)
#        cmd += "--method %s" % method
#        cmd += "-d %s " % txt_db1
        # I would like just raw counts reported
#       if DO_LEN_NORM:
#           cmd += "--lengthNorm "
#       if options.sqlite_db_dir:
#           cmd += "--sqlite_db_dir %s" % options.sqlite_db_dir
#       else: # Using MySQL database
#           if options.passwd == "":
#               cmd += "--host %s --user %s" % (options.host,
#                                               options.user)
#           else:
#               cmd += "--host %s --user %s --passwd %s" % (options.host,
#                                                           options.user,
#                                                           options.passwd)
        os.system(cmd)

        # Add constitutive counts to IR events
        if options.prefix:
            ir_file_name = "%s_all_AS_event_info_irOnly.txt" % prefix
        else:
            ir_file_name = "all_AS_event_info_irOnly.txt"
       
        # IR Counts are already normalized 
        if printExonCoords:
            if os.path.exists(ir_file_name):
                updateCounts2all_as_events(ir_file_name, 
                                       mapped_file1_counts, mapped_file2_counts,
                                       norm1, norm2, jcn_seq_len)

        # The exclusion counts in the original files are incorrect because
        # each end of the junction was calculated separately
        if os.path.exists(ir_file_name):
            fixIRExclusion_count(ir_file_name, all_jcn_count_dict,
                                 norm1, norm2, jcn_seq_len)

        # Combine irOnly with the rest of all_AS_event_info
        if ie1_file is not None:
            if options.prefix:
                main_file = "%s_all_AS_event_info.txt" % prefix
                ir_file = "%s_all_AS_event_info_irOnly.txt" % prefix
            else:
                main_file = "all_AS_event_info.txt"
                ir_file = "all_AS_event_info_irOnly.txt"

            if os.path.exists(ir_file):
                cmd = "cat %s >> %s" % (ir_file, main_file)
                os.system(cmd)

    # Sum Totals for inclusion and exclusion isoforms in all_AS_event_info
    # and adjust for paired end counting
    ## Get sum of exclusion/inclusion for cassette, alternative donor and
    # acceptor, mutually exclusive and multi cassette from acutal files.
    sumExclusion_Inclusion_counts(all_event_info_str, 
                                  ce2total_counts,
                                  alt_donor2total_counts,
                                  alt_accept2total_counts,
                                  afe2total_counts,
                                  ale2total_counts,
                                  mxe2total_counts,
                                  mc2total_counts,
                                  printExonCoords,
                                  norm1, norm2, jcn_seq_len)

    ERROR_LOG.close()                                  

    # Remove all temporary files.
    if not keep_interm:
        os.remove(cassette_out_str)
        os.remove(donor_out_str)
        os.remove(accept_out_str)
        os.remove(jcn_only_donor_out_str)
        os.remove(jcn_only_accept_out_str)
        os.remove(afe_out_str)
        os.remove(ale_out_str)
        os.remove(me_out_str)
        os.remove(mc_out_str)
        os.remove(ir_left_out_str)
        os.remove(ir_right_out_str)
        os.remove(mapped_file1_name)
        os.remove(mapped_file2_name)
        os.remove(exon_coord_file_name)
        os.remove(annot_intron_pick_name)
        os.remove(annot_exon_pick_name)

        if os.path.exists(ir_file):
            os.remove(ir_file)
            if options.prefix:
                os.remove("%s_IR_diff_direction.txt" % prefix)
                os.remove("%s_IR_left_no_partner.txt" % prefix)
                os.remove("%s_IR_right_no_partner.txt" % prefix)
            else:
                os.remove("IR_diff_direction.txt")
                os.remove("IR_left_no_partner.txt")
                os.remove("IR_right_no_partner.txt")

    # Create empty file to indicate that the run completed
    if options.prefix:
        finished_file = open("%s_finished.txt" % prefix, "w")
    else:
        finished_file = open("finished.txt", "w")
    finished_file.write("DONE\n")
    finished_file.close()

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def chrCheck(jcn_chrs, annot_dict, annot_type):
    for chr in jcn_chrs:
        if chr not in annot_dict:
            if annot_type:
                ERROR_LOG.write("WARNING: %s does not exist in %s. Please make sure the genome reference of transcript database matches read alignment reference.\n" % (chr,
                                                                                                                                          annot_type))
            # prevents future key errors
            annot_dict[chr] = [] 


def find_AFE_ALE_clusters(events_dictList,
                          next_or_previous,
                          all_confident_exons, 
                          all_confident_exons_start2end,
                          all_confident_exons_end2start,
                          all_jcn_count_dict, 
                          chr, event_jcns):
    """
    Will update the ad_aa_afe_ale_events which is a list of dictionaries, each
    dictionary will have the following keys:
    isSimple
    isAltFirstLast
    all_coord_end2start/start2end - both can be initiated although only one will
    be functional for each event

    For alternative first and last exons each isoform is represented by the
    most distal jcn (longest intron length)

    If there are also alternative first exon events, there will also be the
    following keys:
    jcn_cluster_sum
    jcn2jcn_str
    jcn2exon_str
    novel_jcns - str
    novel_jcn_sum_raw,
    novel_jcn_sum_lenNorm
    """
    # exon_clusters = [[(jcn,exon), (jcn, exon)],]
    (exon_clusters,
     novel_jcns,
     novel_jcn_sum_samp1,
     novel_jcn_sum_samp2) = find_exon_clusters(next_or_previous,
                                               all_confident_exons,
                                               all_confident_exons_start2end,
                                               all_confident_exons_end2start,
                                               all_jcn_count_dict,
                                               event_jcns)
    if len(exon_clusters) > 1:
        # First event is the AFE/ALE event
        events_dictList[0] =  {"isSimple": False,
                           "isAltFirstLast": True,
                           "all_coord_end2start":{chr:{}},
                           "all_coord_start2end":{chr:{}},
                           "novel_jcns": novel_jcns,
                           "novel_jcn_sum_samp1": novel_jcn_sum_samp1,
                           "novel_jcn_sum_samp2": novel_jcn_sum_samp2,
                           "jcn_cluster_sum":{},
                           "jcn2jcn_str":{},
                           "jcn2exon_str":{}}

        # Treat each cluster as one isoform
        # cluster = [(jcn,exon),]
        for cluster in exon_clusters:
            if len(cluster) == 1:
                chr, start, end = convertCoordStr(cluster[0][0])
                if next_or_previous == "P":
                    updateDictOfLists(events_dictList[0]["all_coord_end2start"][chr],
                                      end, start)
                else:
                    updateDictOfLists(events_dictList[0]["all_coord_start2end"][chr],
                                      start, end)
               
                events_dictList[0]["jcn_cluster_sum"][cluster[0][0]] = all_jcn_count_dict[cluster[0][0]]
                events_dictList[0]["jcn2jcn_str"][cluster[0][0]] = cluster[0][0]
                events_dictList[0]["jcn2exon_str"][cluster[0][0]] = cluster[0][1]
            else:
                # Find longest jcn, 
                # get only exon region that is shared by all exons in cluster
                # get all_jcn_strs and sum
                longest_jcn_len = -1
                longest_jcn = None

                right_most_start = -1
                left_most_end = INFINITY

                jcn_sum = [0,0]
                jcn_strs = []
                for (jcn, exon) in cluster:
                    jcn_strs.append(jcn)
                    jcn_sum[0] += all_jcn_count_dict[jcn][0]
                    jcn_sum[1] += all_jcn_count_dict[jcn][1]

                    chr, jcn_start, jcn_end = convertCoordStr(jcn)
                    this_jcn_len = jcn_end - jcn_start + 1
                    if this_jcn_len > longest_jcn_len:
                        longest_jcn_len = this_jcn_len
                        longest_jcn = jcn

                    if exon:
                        chr, exon_start, exon_end = convertCoordStr(exon)
                        if exon_start > right_most_start:
                            right_most_start = exon_start
                        if exon_end < left_most_end:
                            left_most_end = exon_end

                # Update coord dictionary using shortest
                chr, distal_jcn_start, distal_jcn_end = convertCoordStr(longest_jcn)
                if next_or_previous == "P":
                    updateDictOfLists(events_dictList[0]["all_coord_end2start"][chr],
                                      distal_jcn_end, distal_jcn_start)
                else:
                    updateDictOfLists(events_dictList[0]["all_coord_start2end"][chr],
                                      distal_jcn_start, distal_jcn_end)
                
                events_dictList[0]["jcn_cluster_sum"][longest_jcn] = jcn_sum
                # Making sure that first jcn is prox jcn
                jcn_strs.pop(jcn_strs.index(longest_jcn))
                jcn_strs.insert(0,longest_jcn)
                events_dictList[0]["jcn2jcn_str"][longest_jcn] = ",".join(jcn_strs)
               
                # No overlapping exon found for cluster 
                if right_most_start == -1 or left_most_end < right_most_start:
                    events_dictList[0]["jcn2exon_str"][longest_jcn] = "None"
                else:
                    events_dictList[0]["jcn2exon_str"][longest_jcn] = formatCoordStr(chr, 
                                                                                     right_most_start,
                                                                                     left_most_end)
    else: # Need to remove existing dictionary at first position if there are
          # no more than 1 clusters
        events_dictList.pop(0)
        
    # Now that AFE/ALE event is created, make alt5' and 3' events based on
    # clusters
    # cluster = [(jcn,exon),]
    for cluster in exon_clusters:
        if len(cluster) > 1:
            new_dict = {"isSimple":True,
                        "isAltFirstLast":False,
                        "all_coord_end2start":{chr:{}},
                        "all_coord_start2end":{chr:{}}}

            for jcn, exon in cluster:
                chr, jcn_start, jcn_end = convertCoordStr(jcn)
                if next_or_previous == "P":
                    updateDictOfLists(new_dict["all_coord_end2start"][chr],
                                      jcn_end, jcn_start)
                else: 
                    updateDictOfLists(new_dict["all_coord_start2end"][chr],
                                      jcn_start, jcn_end)

            events_dictList.append(new_dict)

    # No return value

def find_exon_clusters(next_or_previous, all_confident_exons,
                       all_confident_exons_start2end,
                       all_confident_exons_end2start,
                       all_jcn_count_dict, event_jcns):
    """
    (exon_clusters,
     novel_jcns,
     novel_jcn_sum_raw,
     novel_jcn_sum_lenNorm)
    """
    adjacentConfExons = []

    novel_jcns = []
    novel_jcn_sum_raw = 0
    novel_jcn_sum_lenNorm = 0
    # Find all adjacent confident exons
    for jcn_str in event_jcns:
        chr, start, end = convertCoordStr(jcn_str)
        if next_or_previous == "N":
            # chr_start_end
            adjExon = hasAdjExons(chr, all_confident_exons_start2end, [end], "N")
        else:
            adjExon = hasAdjExons(chr, all_confident_exons_end2start, [start], "P")
        if adjExon:
            adjacentConfExons.append((jcn_str, adjExon))
        else: # this is a novel junction
            novel_jcns.append(jcn_str)
            novel_jcn_sum_raw += all_jcn_count_dict[jcn_str][0] 
            novel_jcn_sum_lenNorm += all_jcn_count_dict[jcn_str][1] 
         

    clusters = [[adjacentConfExons.pop()]]

    while len(adjacentConfExons) > 0:
        (this_jcn_str, this_exon) = adjacentConfExons.pop()

        chr, this_exon_start, this_exon_end = convertCoordStr(this_exon)

        overlapping_indices = []

        for i in range(len(clusters)):
            for (c_jcn_str, c_exon) in clusters[i]:
                c_chr, c_exon_start, c_exon_end = convertCoordStr(c_exon)
                if c_chr != chr:
                    print "Error: chromosomes don't match"
                    sys.exit(1)

                if coordsOverlap(this_exon_start, this_exon_end,
                                 c_exon_start, c_exon_end):
                    overlapping_indices.append(i)
                    break

        if len(overlapping_indices) > 1:
            new_cluster = []
            for i in overlapping_indices:
                new_cluster.extend(list(clusters[i]))
            # Add new exon to new_cluster
            new_cluster.append((this_jcn_str, this_exon))

            new_clusters = []
            # Add non overlapping clusters
            for i in range(len(clusters)):
                if i not in overlapping_indices:
                    new_clusters.append(clusters[i])
            # Add the new cluster
            new_clusters.append(new_cluster)
            clusters = new_clusters
        elif len(overlapping_indices) == 1:
            clusters[overlapping_indices[0]].append((this_jcn_str, this_exon))    
        else: # start a new cluster
            clusters.append([(this_jcn_str, this_exon)])

    # Clusters have been found
    return (clusters, ",".join(novel_jcns), novel_jcn_sum_raw,
            novel_jcn_sum_lenNorm)

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def updateCounts2AltDonorAccept(file_out_str, 
                                ir_count_dict,
                                mapped_file1_counts,mapped_file2_counts,
                                norm1, norm2, jcn_seq_len,
                                jcnOnly=False):
    """
    There are a lot of updating to total counts, from the ie junctions and
    exonic junctions.  The total counts reported in the all event file will be
    given in this function.
    """
    
    file = open(file_out_str)
    lines = file.readlines()
    file.close()

    file2 = open(file_out_str, "w")

    # Returned dictionary that will be used later for all_event_info file
    # {(incl_str, excl_str): (excl1, incl1, excl2, incl2)}
    event2counts = {}
    

    for line in lines:
        line = formatLine(line)

        line_list = line.split("\t")

        chr = line_list[3]

        excl_raw_str = line_list[-5]
        jcn_incl_raw = int(line_list[-4])
   
        incl_add_coord = line_list[-1]

        # Find all inclusion regions
#        incl_junction = .join([line_list[3], line_list[5], line_list[6]])
        incl_junction = formatCoordStr(line_list[3], line_list[5], line_list[6])

        incl_start = int(line_list[5])
        incl_end = int(line_list[6]) 
        excl_jcns = line_list[7].split(";")

        exclusion_raw_list = map(int, excl_raw_str.split(";"))

        event_key = (incl_junction,
                     ";".join(excl_jcns))


        alt_start_or_end = determineAltStartOrEnd(incl_start, incl_end,
                                                  excl_jcns)
       
        (ordered_pos, 
         proportions1, 
         proportions2,
         total_ordered_raw_list,
         total_ordered_lenNorm_list) = getSSOrderAndProportions(alt_start_or_end,
                                                  incl_start, incl_end,
                                                  excl_jcns,
                                                  exclusion_raw_list,
                                                  jcn_incl_raw,
                                                  exclusion_raw_list,
                                                  jcn_incl_raw)

        # Reset values of total_ordered_lenNorm_list
        for i in range(len(total_ordered_lenNorm_list)):
            total_ordered_lenNorm_list[i] = 0

        if jcnOnly:
            incl_raw = jcn_incl_raw

            # Add the rest of the exclusion junction counts to the total exclusion
            excl_raw = sum(exclusion_raw_list)

            if hasNegativeVals(excl_raw, incl_raw, excl_raw, incl_raw):
                ERROR_LOG.write("Negative Vals: %s\n" % line)
                excl_raw = 0
                incl_raw = 0
                excl_lenNorm = 0
                incl_lenNorm = 0

            out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                                excl_raw, incl_raw,
                                                excl_raw, incl_raw)

            file2.write(out_str)
            continue



        # Find IE Junction counts
        ie_jcn_cts_raw = []

        if alt_start_or_end == "alt_start":
            for i in range(len(ordered_pos)-1):
                intron = formatCoordStr(chr, ordered_pos[i], incl_end)
           
                ie_jcn = formatCoordStr(chr, ordered_pos[i] - 1,
                                       ordered_pos[i]) 

                ie_jcn_ct_raw = 0

                if ir_count_dict:
                    if intron in ir_count_dict:

                        ie_jcn_ct_raw = ir_count_dict[intron]["left"][1]

                        if norm2:
                            ie_jcn_ct_raw = int(round(ie_jcn_ct_raw/norm2))

                ie_jcn_cts_raw.append(ie_jcn_ct_raw)

        else: # alt_end
            for i in range(1,len(ordered_pos)):
                intron = formatCoordStr(chr, incl_start, ordered_pos[i])

                ie_jcn = formatCoordStr(chr, ordered_pos[i],
                                       ordered_pos[i] + 1)

                ie_jcn_ct_raw = 0

                if ir_count_dict:
                    if intron in ir_count_dict:
                        ie_jcn_ct_raw = ir_count_dict[intron]["right"][1]

                        if norm2:
                            ie_jcn_ct_raw = int(round(ie_jcn_ct_raw/norm2))

                ie_jcn_cts_raw.append(ie_jcn_ct_raw)

        ### Split IE junction counts by each isoform
        (sub_proportions1,
         sub_proportions2) = updateProportions(alt_start_or_end,
                                               proportions1, proportions2)
        
        (ie_isoform_cts1,
         ie_isoform_cts2) =  splitSharedRegions(alt_start_or_end,
                                                ie_jcn_cts_raw,
                                                ie_jcn_cts_raw,
                                                sub_proportions1,
                                                sub_proportions2)

        # parallel to ordered_pos
        isoform_lengths = getAD_AA_isoform_lengths(alt_start_or_end, 
                                                   ordered_pos, jcn_seq_len)

        # Add counts
        incl_raw = 0
        incl_lenNorm = 0
        excl_raw = 0
        excl_lenNorm = 0
        
        incl_ordered_pos = None
        if alt_start_or_end == "alt_start":
            incl_ordered_pos_idx = ordered_pos.index(incl_start)
        else:
            incl_ordered_pos_idx = ordered_pos.index(incl_end)

        # Normalize existing junction counts and add to total
        for i in range(len(ordered_pos)):
            update_ct_raw = total_ordered_raw_list[i]
            total_ordered_raw_list[i] = update_ct_raw
            update_ct_lenNorm = normalizeByLen(total_ordered_raw_list[i], 
                                        isoform_lengths[i])
            total_ordered_lenNorm_list[i] = update_ct_lenNorm

            if i == incl_ordered_pos_idx:
                incl_raw += update_ct_raw
                incl_lenNorm += update_ct_lenNorm
            else:
                excl_raw += update_ct_raw
                excl_lenNorm += update_ct_lenNorm


        if alt_start_or_end == "alt_start":
            # When the longest intron is the "inclusion" isoform, all other
            # counts are added to the exclusion counts
            if ordered_pos.index(incl_start) == 0:
                for i in range(len(ie_isoform_cts2)):
                    this_ie_ct_raw = ie_isoform_cts2[i]
                    this_ie_ct_lenNorm = normalizeByLen(ie_isoform_cts2[i],
                                                 isoform_lengths[i+1])

                    total_ordered_raw_list[i+1] += this_ie_ct_raw
                    total_ordered_lenNorm_list[i+1] += this_ie_ct_lenNorm

                    excl_raw += this_ie_ct_raw
                    excl_lenNorm += this_ie_ct_lenNorm
                                                                   
            else:
                for i in range(len(ie_isoform_cts2)):
                    this_ie_ct_raw = ie_isoform_cts2[i]
                    this_ie_ct_lenNorm = normalizeByLen(ie_isoform_cts2[i],
                                                 isoform_lengths[i+1])
                
                    total_ordered_raw_list[i+1] += this_ie_ct_raw
                    total_ordered_lenNorm_list[i+1] += this_ie_ct_lenNorm

                    if i == (ordered_pos.index(incl_start) - 1):
                        incl_raw += this_ie_ct_raw
                        incl_lenNorm += this_ie_ct_lenNorm
                    else:
                        excl_raw += this_ie_ct_raw
                        excl_lenNorm += this_ie_ct_lenNorm
          
        else: # "alt_end"
            # When the longest intron is the "inclusion" isoform, all other
            # counts are added to the exclusion counts
            if ordered_pos.index(incl_end) == (len(ordered_pos) - 1):
                for i in range(len(ie_isoform_cts2)):
                    this_ie_ct_raw = ie_isoform_cts2[i]
                    this_ie_ct_lenNorm = normalizeByLen(ie_isoform_cts2[i],
                                                 isoform_lengths[i])

                    total_ordered_raw_list[i] += this_ie_ct_raw
                    total_ordered_lenNorm_list[i] += this_ie_ct_lenNorm
                    
                    excl_raw += this_ie_ct_raw
                    excl_lenNorm += this_ie_ct_lenNorm

            else:
                for i in range(len(ie_isoform_cts2)):
                    this_ie_ct_raw = ie_isoform_cts2[i]
                    this_ie_ct_lenNorm = normalizeByLen(ie_isoform_cts2[i],
                                                 isoform_lengths[i])

                    total_ordered_raw_list[i] += this_ie_ct_raw
                    total_ordered_lenNorm_list[i] += this_ie_ct_lenNorm

                    if i == ordered_pos.index(incl_end):
                        incl_raw += this_ie_ct_raw
                        incl_lenNorm += this_ie_ct_lenNorm
                    else:
                        excl_raw += this_ie_ct_raw
                        excl_lenNorm += this_ie_ct_lenNorm


        #### Now split up the inclusion add coords by the proportions.
        incl_regions = incl_add_coord.split(";")

        mapped_file1_incl_counts = []
        mapped_file2_incl_counts = []

        for incl_region in incl_regions:
            if incl_region in mapped_file1_counts:
                mapped_file1_incl_counts.append(mapped_file1_counts[incl_region])  
            else:
                mapped_file1_incl_counts.append(0)

            if incl_region in mapped_file2_counts:
                mapped_file2_incl_counts.append(mapped_file2_counts[incl_region])
            else:
                mapped_file2_incl_counts.append(0)
    
        # Inclusion regions are added proportionally based on junction
        # proportions
        # The inclusion regions are shared by a subset of the isoforms,
        # not all; therefore, the proportions need to be recalculated
        (exonic_isoform_cts1,
         exonic_isoform_cts2) = splitSharedRegions(alt_start_or_end,
                                            mapped_file1_incl_counts,
                                            mapped_file2_incl_counts,
                                            sub_proportions1, 
                                            sub_proportions2)

        if alt_start_or_end == "alt_start":
            # If the inclusion isoform is the longest intron, then all counts
            # go to the exclusion isoform
            if ordered_pos.index(incl_start) == 0:
                for i in range(len(exonic_isoform_cts2)):
                    this_exon_ct_raw = exonic_isoform_cts2[i]
                    this_exon_ct_lenNorm = normalizeByLen(exonic_isoform_cts2[i],
                                                   isoform_lengths[i+1])

                    total_ordered_raw_list[i+1] += this_exon_ct_raw
                    total_ordered_lenNorm_list[i+1] += this_exon_ct_lenNorm
                    
                    excl_raw += this_exon_ct_raw
                    excl_lenNorm += this_exon_ct_lenNorm
            else:
                for i in range(len(exonic_isoform_cts2)):
                    this_exon_ct_raw = exonic_isoform_cts2[i]
                    this_exon_ct_lenNorm = normalizeByLen(exonic_isoform_cts2[i],
                                                   isoform_lengths[i+1])

                    total_ordered_raw_list[i+1] += this_exon_ct_raw
                    total_ordered_lenNorm_list[i+1] += this_exon_ct_lenNorm

                    if i == (ordered_pos.index(incl_start) - 1):
                        incl_raw += this_exon_ct_raw
                        incl_lenNorm += this_exon_ct_lenNorm
                    else:
                        excl_raw += this_exon_ct_raw
                        excl_lenNorm += this_exon_ct_lenNorm

        else: # alt_end
            if ordered_pos.index(incl_end) == (len(ordered_pos) - 1):
                for i in range(len(exonic_isoform_cts2)):
                    this_exon_ct_raw = exonic_isoform_cts2[i]
                    this_exon_ct_lenNorm = normalizeByLen(exonic_isoform_cts2[i],
                                                    isoform_lengths[i])

                    total_ordered_raw_list[i] += this_exon_ct_raw
                    total_ordered_lenNorm_list[i] += this_exon_ct_lenNorm

                    excl_raw += this_exon_ct_raw
                    excl_lenNorm += this_exon_ct_lenNorm
            else:
                for i in range(len(exonic_isoform_cts2)):
                    this_exon_ct_raw = exonic_isoform_cts2[i]
                    this_exon_ct_lenNorm = normalizeByLen(exonic_isoform_cts2[i],
                                                   isoform_lengths[i])            

                    total_ordered_raw_list[i] += this_exon_ct_raw
                    total_ordered_lenNorm_list[i] += this_exon_ct_lenNorm

                    if i == ordered_pos.index(incl_end):
                        incl_raw += this_exon_ct_raw
                        incl_lenNorm += this_exon_ct_lenNorm
                    else:
                        excl_raw += this_exon_ct_raw
                        excl_lenNorm += this_exon_ct_lenNorm
    
        # Create exclusion cts list
        if alt_start_or_end == "alt_start":
            this_incl_raw_ct = total_ordered_raw_list.pop(ordered_pos.index(incl_start))
            this_incl_lenNorm_ct = total_ordered_lenNorm_list.pop(ordered_pos.index(incl_start))
        else: # alt_end
            this_incl_raw_ct = total_ordered_raw_list.pop(ordered_pos.index(incl_end))
            this_incl_lenNorm_ct = total_ordered_lenNorm_list.pop(ordered_pos.index(incl_end))

        if this_incl_raw_ct != incl_raw:
            ERROR_LOG.write("updateCounts2AltDonorAccept: inclusion raw counts do not agree: %s\n" % event_key)
        if this_incl_lenNorm_ct != incl_lenNorm:
            ERROR_LOG.write("updateCounts2AltDonorAccept: inclusion lenNorm counts do not agree: %s\n" % event_key)

        total_excl_raw_list = total_ordered_raw_list
        total_excl_lenNorm_list = total_ordered_lenNorm_list

#       (ordered_pos, 
#        new_proportions1, 
#        new_proportions2,
#        not_used1,
#        not_used2) = getSSOrderAndProportions(alt_start_or_end,
#                                              incl_start, incl_end,
#                                              excl_jcns,
#                                              total_excl_cts1_list,
#                                              incl1,
#                                              total_excl_cts2_list,
#                                              incl2)

#       e_or_i = checkExclusionInclusion_AA_AD_AFE_ALE(alt_start_or_end,
#                                                      incl_start, incl_end,
#                                                      ordered_pos,
#                                                      new_proportions1,
#                                                      new_proportions2)
#       line_list[1] = e_or_i

        if hasNegativeVals(excl_raw, incl_raw, excl_lenNorm, incl_lenNorm):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            excl_raw = 0
            incl_raw = 0
            excl_lenNorm = 0
            incl_lenNorm = 0

        out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                            excl_raw, incl_raw,
                                            excl_lenNorm, incl_lenNorm)

        
        # Add counts to event dictionary
        event2counts[event_key] = (excl_raw, incl_raw, excl_lenNorm, incl_lenNorm)

        file2.write(out_str)

    file2.close()

    return event2counts

def getAD_AA_isoform_lengths(alt_start_or_end, ordered_pos, jcn_seq_len):

    isoform_lengths = []
    if alt_start_or_end == "alt_start":
        # First pos, is just junction
        isoform_lengths.append(jcn_seq_len)
        for i in range(len(ordered_pos)-1):
            this_exonic_portion = ordered_pos[i+1] - ordered_pos[0]
            num_ie_jcns = i + 1
            # isoform length is exonic portion + jcn_seq_len + (num_ie_jcns)*ie_jcn_len
            isoform_lengths.append(this_exonic_portion + jcn_seq_len + (num_ie_jcns*jcn_seq_len))
    else: # alt_end
        for i in range(len(ordered_pos)-1):
            this_exonic_portion = ordered_pos[-1] - ordered_pos[i]
            num_ie_jcns = len(ordered_pos) - i -  1
            isoform_lengths.append(this_exonic_portion + jcn_seq_len + (num_ie_jcns*jcn_seq_len))
        # Last position is just jcn sequence
        isoform_lengths.append(jcn_seq_len)

    return isoform_lengths
            

def updateCounts2AFE_ALE(a_out_str, 
                         all_jcn_count_dict,
                         all_coord_start2end,
                         all_coord_end2start,
                         mapped_file1_counts, mapped_file2_counts,
                         norm1, norm2, jcn_seq_len):

    """
    For length normalization:
    - Make a dictionary associating distal junction to isoform length
    - build it by associating distal junctions to exonic coords
    """

    file = open(a_out_str)
    lines = file.readlines()
    file.close()

    file2 = open(a_out_str, "w")

    # Returned dictionary that will be used later for all_event_info file
    # {(incl_str, excl_str): (excl1, incl1, excl2, incl2)}
    event2counts = {}

    for line in lines:
        line = formatLine(line)

        line_list = line.split("\t")

        # Used for length normalization
        distal_jcn2totalCounts_raw = {}
        distal_jcn2totalCounts_lenNorm = {}
        distal_jcn2isoform_length = {}

        distal_jcn2excl_ct_list_idx = {}

        incl_add_coord = line_list[11]

        incl_raw = 0
        incl_lenNorm = 0

        excl_raw_ct_list = map(int, line_list[7].split(";"))
        excl_lenNorm_ct_list = map(int, line_list[9].split(";"))

        excl_raw = 0
        excl_lenNorm = 0

        incl_jcns = line_list[5]
        excl_jcns_full = line_list[6].split(";")
        excl_jcns = []

        for i in range(len(excl_jcns_full)):
            # First exclusion junction is the distal junction, representing
            # the group
            jcns_in_group = excl_jcns_full[i].split(",")
            this_distal = jcns_in_group[0]

            excl_jcns.append(this_distal)
            distal_jcn2totalCounts_raw[this_distal] = excl_raw_ct_list[i]
            distal_jcn2totalCounts_lenNorm[this_distal] = excl_lenNorm_ct_list[i]
            # Working out the math, I only normalize by one junction sequence
            # length
            distal_jcn2isoform_length[this_distal] = jcn_seq_len
            distal_jcn2excl_ct_list_idx[this_distal] = i

        event_key = (incl_jcns,
                     line_list[6])             

        incl_jcns_in_group = incl_jcns.split(",")
        distal_jcn = incl_jcns_in_group[0]
        distal_jcn2totalCounts_raw[distal_jcn] = int(line_list[8])
        distal_jcn2totalCounts_lenNorm[distal_jcn] = int(line_list[10])
        distal_jcn2isoform_length[distal_jcn] = jcn_seq_len

        chr, incl_start, incl_end = convertCoordStr(distal_jcn)

        alt_start_or_end = determineAltStartOrEnd(incl_start, incl_end,
                                                  excl_jcns)

        # Inclusion counts
        prop1 = None
        prop2 = None
        if alt_start_or_end == "alt_start":
            (prop1,
             prop2) = getJunctionProportion_AFE_ALE(chr,
                                           incl_jcns,
                                           all_jcn_count_dict,
                                           "start",
                                           all_coord_start2end)
        else: # alt_end
            (prop1,
             prop2) = getJunctionProportion_AFE_ALE(chr,
                                           incl_jcns,
                                           all_jcn_count_dict,
                                           "end",
                                           all_coord_end2start)

        if incl_add_coord != "None":
            this_chr, incl_exon_start, incl_exon_end = convertCoordStr(incl_add_coord)
            incl_exon_len = incl_exon_end - incl_exon_start + 1
            this_distal_jcn = findDistalJcn(distal_jcn2isoform_length.keys(), 
                                            incl_add_coord)
            if this_distal_jcn:
                distal_jcn2isoform_length[this_distal_jcn] += incl_exon_len
            else:
                ERROR_LOG.write("Can't find associated junction for AFE/ALE exon %s\n" % incl_add_coord)

            if incl_add_coord in mapped_file2_counts:
                if this_distal_jcn:
                    distal_jcn2totalCounts_raw[this_distal_jcn] += int(round(mapped_file2_counts[incl_add_coord] * prop2))

#           if incl_add_coord in mapped_file2_counts:
#               if this_distal_jcn:
#                   distal_jcn2totalCounts2[this_distal_jcn] += int(round(mapped_file2_counts[incl_add_coord] * prop2))

        # Exclusion counts to all other exons
        excl_add_coords = line_list[12]
        excl_add_coord_list = []
        total_excl_raw_list = list(excl_raw_ct_list)
        total_excl_lenNorm_list = list(excl_lenNorm_ct_list)


        excl_raw_add = 0
#        excl2_add = 0
        if excl_add_coords != "None":
            excl_add_coord_list = excl_add_coords.split(",")

            for excl_add_coord in excl_add_coord_list:
                this_distal_jcn = findDistalJcn(distal_jcn2isoform_length.keys(), 
                                                excl_add_coord)

                if not this_distal_jcn:
                    ERROR_LOG.write("Can't find associated junction for AFE/ALE exon %s\n" % excl_add_coord)
                    continue

                chr, excl_start, excl_end = convertCoordStr(this_distal_jcn)
                chr, this_c_start, this_c_end = convertCoordStr(excl_add_coord)
    
                distal_jcn2isoform_length[this_distal_jcn] += (this_c_end - this_c_start + 1)

                if alt_start_or_end == "alt_start":
                    if excl_start - 1 != this_c_end:
                        continue
                else:
                    if excl_end + 1 != this_c_start:
                        continue

                # Find which set of exclusion junctions this distal
                # exclusion junction is from
                excl_jcn_set = findExclSet(this_distal_jcn, excl_jcns_full)
    
                prop1 = None
                prop2 = None
                if alt_start_or_end == "alt_start":
                    (prop1,
                     prop2) = getJunctionProportion_AFE_ALE(chr,
                                                   excl_jcn_set,
                                                   all_jcn_count_dict,
                                                   "start",
                                                   all_coord_start2end)
                else: # alt_end
                    (prop1,
                     prop2) = getJunctionProportion_AFE_ALE(chr,
                                                   excl_jcn_set,
                                                   all_jcn_count_dict,
                                                   "end",
                                                   all_coord_end2start)

                this_idx = excl_jcns.index(this_distal_jcn)

#               if excl_add_coord in mapped_file1_counts:
#                   distal_jcn2totalCounts1[this_distal_jcn] += int(round(mapped_file1_counts[excl_add_coord] * prop1))
#                   total_excl_cts1_list[this_idx] += int(round(mapped_file1_counts[excl_add_coord] * prop1))
                if excl_add_coord in mapped_file2_counts:
                    distal_jcn2totalCounts_raw[this_distal_jcn] += int(round(mapped_file2_counts[excl_add_coord] * prop2))
                    total_excl_raw_list[this_idx] += int(round(mapped_file2_counts[excl_add_coord] * prop2))

        # Length normalize
        incl_raw = distal_jcn2totalCounts_raw[distal_jcn]
        incl_lenNorm = normalizeByLen(distal_jcn2totalCounts_raw[distal_jcn],
                                      distal_jcn2isoform_length[distal_jcn])


        for excl_jcn in distal_jcn2excl_ct_list_idx:
            excl_raw += distal_jcn2totalCounts_raw[excl_jcn]
            excl_lenNorm += normalizeByLen(distal_jcn2totalCounts_raw[excl_jcn],
                                    distal_jcn2isoform_length[excl_jcn])

            total_excl_raw_list[distal_jcn2excl_ct_list_idx[excl_jcn]] = distal_jcn2totalCounts_raw[excl_jcn]
            total_excl_lenNorm_list[distal_jcn2excl_ct_list_idx[excl_jcn]] = normalizeByLen(distal_jcn2totalCounts_raw[excl_jcn],
                                                                                            distal_jcn2isoform_length[excl_jcn])


#       (ordered_pos, 
#        new_proportions1, 
#        new_proportions2,
#        not_used1,
#        not_used2) = getSSOrderAndProportions(alt_start_or_end,
#                                              incl_start, incl_end,
#                                              excl_jcns,
#                                              total_excl_cts1_list,
#                                              incl1,
#                                              total_excl_cts2_list,
#                                              incl2)

#       e_or_i = checkExclusionInclusion_AA_AD_AFE_ALE(alt_start_or_end,
#                                                      incl_start, incl_end,
#                                                      ordered_pos,
#                                                      new_proportions1,
#                                                      new_proportions2)

#       line_list[1] = e_or_i

        if hasNegativeVals(excl_raw, incl_raw, excl_lenNorm, incl_lenNorm):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            excl_raw = 0
            incl_raw = 0
            excl_lenNorm = 0
            incl_lenNorm = 0

        out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                            excl_raw, incl_raw, 
                                            excl_lenNorm, incl_lenNorm)
        # Add counts to even dictionary
        event2counts[event_key] = (excl_raw, incl_raw, excl_lenNorm, incl_lenNorm)

        file2.write(out_str)

    file2.close()

    return event2counts

def findDistalJcn(jcn_list, incl_add_coord):
    chr, exon_start, exon_end = convertCoordStr(incl_add_coord)

    for jcn in jcn_list:
        chr, intron_start, intron_end = convertCoordStr(jcn)

        if (intron_start - 1) == exon_end:
            return jcn
        if (intron_end + 1) == exon_start:
            return jcn

    # If not returned, could not find adjacent exon
    return None


def updateCounts2Cassette(file_out_str, 
                          mapped_file1_counts, mapped_file2_counts,
                          norm1, norm2, jcn_seq_len):
    
    file = open(file_out_str)
    lines = file.readlines()
    file.close()

    file2 = open(file_out_str, "w")

    # {exon_str: (excl1, incl1, excl2, incl2)}
    ce2total_counts = {}

    for line in lines:
        line = formatLine(line)

        line_list = line.split("\t")

        upstrm_jcns = line_list[5].split(",")
        dwnstrm_jcns = line_list[8].split(",")

        excl_jcns = inferExclusionJunctions(upstrm_jcns, dwnstrm_jcns)

        excl_raw = int(line_list[-5])
        incl_raw = int(line_list[-4])
        excl_lenNorm = int(line_list[-3])
        incl_lenNorm = int(line_list[-2])
    
        incl_add_coord = line_list[-1]

        # Normalize inclusion counts by exon length
        ex_chr, ex_start, ex_end = convertCoordStr(incl_add_coord)
        exon_len = ex_end - ex_start + 1
        incl_len = (2 * jcn_seq_len) + exon_len

        if incl_add_coord in mapped_file2_counts:
            incl_raw += mapped_file2_counts[incl_add_coord]
            incl_lenNorm += normalizeByLen(mapped_file2_counts[incl_add_coord], incl_len)

#       e_or_i = checkExclusionInclusion(excl1,
#                                       incl1,
#                                       excl2,
#                                       incl2) 

        if hasNegativeVals(excl_raw, incl_raw, excl_lenNorm, incl_lenNorm):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            excl_raw = 0
            incl_raw = 0
            excl_lenNorm = 0
            incl_lenNorm = 0

        out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                            excl_raw, incl_raw,
                                            excl_lenNorm, incl_lenNorm)

        ce2total_counts[incl_add_coord] = (excl_raw, incl_raw, excl_lenNorm, incl_lenNorm)

        file2.write(out_str)

    file2.close()

    return ce2total_counts 

def updateCounts2all_as_events(file_str, 
                               mapped_file1_counts, mapped_file2_counts,
                               norm1, norm2, jcn_seq_len):

    """
    Update exon counts
    """

    file = open(file_str)
    lines = file.readlines()
    file.close()

    file2 = open(file_str, "w")

    for line in lines:
        line = line.rstrip("\n")
    
        line_list = line.split("\t")

        excl_exons = line_list[7]
        incl_exons = line_list[8]

        const_exons = line_list[10]

        (excl_ct_str_samp1, excl_ct_str_samp2,
         sum_excl_ct_samp1, sum_excl_ct_samp2) = getCoordCounts4all_as_events(excl_exons, 
                                                                              mapped_file1_counts,
                                                                              mapped_file2_counts)
        (incl_ct_str_samp1, incl_ct_str_samp2,
         sum_incl_ct_samp1, sum_incl_ct_samp2) = getCoordCounts4all_as_events(incl_exons, 
                                                                              mapped_file1_counts,
                                                                              mapped_file2_counts)

        (const_ct_str_samp1, const_ct_str_samp2,
         sum_const_ct_samp1, sum_const_ct_samp2) = getCoordCounts4all_as_events(const_exons, 
                                                                                mapped_file1_counts,
                                                                                mapped_file2_counts)

#       excl_jcns = parse_all_as_event_regions(line_list[5])
#       incl_jcns = parse_all_as_event_regions(line_list[6])

#       excl_exons = parse_all_as_event_regions(line_list[7])
#       incl_exons = parse_all_as_event_regions(line_list[8])

#       ie_jcns = parse_all_as_event_regions(line_list[9])

        if sum_excl_ct_samp1 is None:   
            sum_excl_ct_samp1 = 0
            sum_excl_ct_samp2 = 0

        if sum_incl_ct_samp1 is None:
            sum_incl_ct_samp1 = 0
            sum_incl_ct_samp2 = 0

        if sum_const_ct_samp1 is None:
            sum_const_ct_samp1 = 0
            sum_const_ct_samp2 = 0

        # Length normalize exonic counts constitutive regions are not used in
        # length normalization
#       # Add jcns first
#       incl_isoform_length = jcn_seq_len * (len(incl_jcns) + len(ie_jcns))
#       excl_isoform_length = jcn_seq_len * len(excl_jcns)

#       # Add exon lengths
#       for exon_coord in incl_exons:
#           if exon_coord != "" and exon_coord != "None":
#               this_chr, exon_start, exon_end = convertCoordStr(exon_coord)
#               incl_isoform_length += (exon_end - exon_start + 1)
#       for exon_coord in excl_exons:
#           if exon_coord != "" and exon_coord != "None":
#               this_chr, exon_start, exon_end = convertCoordStr(exon_coord)
#               excl_isoform_length += (exon_end - exon_start + 1)

#       sum_excl_ct_raw = normalizeByLen(sum_excl_ct_raw,
#                                          excl_isoform_length)
#       sum_excl_ct_lenNorm = normalizeByLen(sum_excl_ct_lenNorm,
#                                          excl_isoform_length)
#       sum_incl_ct_raw = normalizeByLen(sum_incl_ct_raw,
#                                          incl_isoform_length)
#       sum_incl_ct_lenNorm = normalizeByLen(sum_incl_ct_lenNorm,
#                                          incl_isoform_length)

        out_str = getAllEventStr(line_list[0], 
                                 line_list[1],
                                 line_list[2],
                                 line_list[3],
                                 line_list[4],
                                 line_list[5],
                                 line_list[6],
                                 line_list[7],
                                 line_list[8],
                                 line_list[9],
                                 line_list[10],
                                 line_list[11],
                                 line_list[12],
                                 line_list[13],
                                 line_list[14],
                                 excl_ct_str_samp2,
                                 incl_ct_str_samp2,
                                 sum_excl_ct_samp2,
                                 sum_incl_ct_samp2,
                                 line_list[19],
                                 line_list[20],
                                 const_ct_str_samp2,
                                 sum_const_ct_samp2)

        file2.write(out_str + "\n")

    file2.close()

def updateCounts2MutuallyExclusive(me_out_str, 
                                   all_jcn_count_dict,
                                   all_coord_start2end,
                                   all_coord_end2start,
                                   mapped_file1_counts, mapped_file2_counts,
                                   norm1, norm2, jcn_seq_len):
    file = open(me_out_str)
    lines = file.readlines()
    file.close()

    file2 = open(me_out_str, "w")

    # returned dictionary. Holds the total counts.
    # {(incl_exon, excl_exons): (excl1, incl1, excl2, incl2)}
    event2counts = {}

    for line in lines:
        line = formatLine(line)

        line_list = line.split("\t")

        excl_raw = int(line_list[9])
        incl_raw = int(line_list[10])
        excl_lenNorm = int(line_list[11])
        incl_lenNorm = int(line_list[12])

        incl_add_coord = line_list[13]

        chr, incl_exon_start, incl_exon_end = convertCoordStr(incl_add_coord)
        incl_isoform_len = (incl_exon_end - incl_exon_start + 1) + (2*jcn_seq_len)

        upstrm_jcn_start = int(line_list[5])
        dwnstrm_jcn_end = int(line_list[6])

        first_incl_jcn = formatCoordStr(chr, upstrm_jcn_start, 
                                            incl_exon_start - 1)
        second_incl_jcn = formatCoordStr(chr, incl_exon_end + 1,
                                        dwnstrm_jcn_end)
        incl_jcns = [first_incl_jcn, second_incl_jcn]
        
        excl_add_coords = line_list[14]

        excl_add_coord_list = excl_add_coords.split(",")

        event_key = (incl_add_coord, excl_add_coords.replace(",",";"))

        excl_jcns = []
        for excl_add_coord in excl_add_coord_list:
            chr, excl_exon_start, excl_exon_end = convertCoordStr(excl_add_coord)
            first_excl_jcn = formatCoordStr(chr, upstrm_jcn_start,
                                           excl_exon_start - 1)
            second_excl_jcn = formatCoordStr(chr, excl_exon_end + 1,
                                            dwnstrm_jcn_end)

            excl_jcns.append(first_excl_jcn)
            excl_jcns.append(second_excl_jcn)

        # Inclusion counts
        (mxe_proportion1,
         mxe_proportion2) = getMXE_MCProportion(incl_add_coord,
                                             upstrm_jcn_start,
                                             dwnstrm_jcn_end,
                                             all_jcn_count_dict,
                                             all_coord_start2end,
                                             all_coord_end2start)

#       if incl_add_coord in mapped_file1_counts:
#           incl1 += normalizeByLen(int(round(mapped_file1_counts[incl_add_coord] * mxe_proportion1)),
#                                   incl_isoform_len)
        if incl_add_coord in mapped_file2_counts:
            incl_raw += int(round(mapped_file2_counts[incl_add_coord] * mxe_proportion2))
            incl_lenNorm += normalizeByLen(int(round(mapped_file2_counts[incl_add_coord] * mxe_proportion2)),
                                    incl_isoform_len)

        # Exclusion counts to all other exons
        excl_raw_add = 0
        excl_lenNorm_add = 0
        for excl_add_coord in excl_add_coord_list:

            this_chr, excl_exon_start, excl_exon_end = convertCoordStr(excl_add_coord)
            this_isoform_len = (excl_exon_end - excl_exon_start + 1) + (2*jcn_seq_len)

            (mxe_proportion1,
             mxe_proportion2) = getMXE_MCProportion(excl_add_coord,
                                                 upstrm_jcn_start,
                                                 dwnstrm_jcn_end,
                                                 all_jcn_count_dict,
                                                 all_coord_start2end,
                                                 all_coord_end2start)

#           if excl_add_coord in mapped_file1_counts:
#               excl1_add += normalizeByLen(int(round(mapped_file1_counts[excl_add_coord]* mxe_proportion1)),
#                                           this_isoform_len)
            if excl_add_coord in mapped_file2_counts:
                excl_raw_add += int(round(mapped_file2_counts[excl_add_coord]* mxe_proportion2))
                excl_lenNorm_add += normalizeByLen(int(round(mapped_file2_counts[excl_add_coord]* mxe_proportion2)),
                                            this_isoform_len)

        excl_raw += excl_raw_add
        excl_lenNorm += excl_lenNorm_add

#       e_or_i = checkExclusionInclusion(excl1,
#                                       incl1,
#                                       excl2,
#                                       incl2) 

#       line_list[1] = e_or_i

        if hasNegativeVals(excl_raw, incl_raw, excl_lenNorm, incl_lenNorm):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            excl_raw = 0
            incl_raw = 0
            excl_lenNorm = 0
            incl_lenNorm = 0

        out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                           excl_raw, incl_raw, 
                                           excl_lenNorm, incl_lenNorm)
        file2.write(out_str)

        event2counts[event_key] = (excl_raw, incl_raw, excl_lenNorm, incl_lenNorm)

    file2.close()
    
    return event2counts

def updateCounts2MultiCassette(mc_out_str, 
                               all_jcn_count_dict,
                               all_coord_start2end,
                               all_coord_end2start,
                               mapped_file1_counts, mapped_file2_counts,
                               norm1, norm2, jcn_seq_len):
    file = open(mc_out_str)
    lines = file.readlines()
    file.close()

    file2 = open(mc_out_str, "w")

    # returned dictionary. Holds the total counts.
    # {(incl_exons, excl_jcn): (excl1, incl1, excl2, incl2)}
    event2counts = {}

    for line in lines:
        line = formatLine(line)

        line_list = line.split("\t")

        excl_raw = int(line_list[8])
        incl_raw = int(line_list[9])
        excl_lenNorm = int(line_list[10])
        incl_lenNorm = int(line_list[11])

        incl_add_coords = line_list[12]
        incl_add_coord_list = incl_add_coords.split(",")

        inclusion_len = 0
        for incl_add_coord in incl_add_coord_list:
            this_chr, this_start, this_end = convertCoordStr(incl_add_coord)
            inclusion_len += (this_end - this_start + 1)

        chr = line_list[3]
        excl_start = int(line_list[5])
        excl_end = int(line_list[6])

        incl_jcns = inferInclusionJunctions(chr, excl_start, excl_end, incl_add_coord_list)
        # Add the number of jcn sequence to inclusion isoform length
        inclusion_len += (len(incl_jcns) * jcn_seq_len)
        
        excl_jcn = formatCoordStr(chr, excl_start, excl_end)

        event_key = (incl_add_coords.replace(",",";"), excl_jcn)

        incl_jcn_coords = []
        for incl_jcn in incl_jcns:
            chr, start, end = convertCoordStr(incl_jcn)
            incl_jcn_coords.append((start,end))

        incl_jcn_coords.sort()
        
        # Loop through every exon and add counts
        last_start = excl_start
        last_end = incl_jcn_coords[0][1]

        incl_raw_add = 0
        incl_lenNorm_add = 0
        for incl_jcn_coord in incl_jcn_coords[1:]:
            this_start = incl_jcn_coord[0]
            this_end = incl_jcn_coord[1]

            inferred_exon = formatCoordStr(chr,
                                          last_end + 1,
                                          this_start - 1)

            (mc_proportion1,
             mc_proportion2) = getMXE_MCProportion(inferred_exon,
                                                   last_start,
                                                   this_end,
                                                   all_jcn_count_dict,
                                                   all_coord_start2end,
                                                   all_coord_end2start)

#           if inferred_exon in mapped_file1_counts:
#               incl1_add += normalizeByLen(int(round(mapped_file1_counts[inferred_exon] * mc_proportion1)),
#                                           inclusion_len)

            if inferred_exon in mapped_file2_counts:
                incl_raw_add += int(round(mapped_file2_counts[inferred_exon] * mc_proportion2))
                incl_lenNorm_add += normalizeByLen(int(round(mapped_file2_counts[inferred_exon] * mc_proportion2)),
                                            inclusion_len)

            last_end = this_end
            last_start = this_start

        incl_raw += incl_raw_add
        incl_lenNorm += incl_lenNorm_add

#       e_or_i = checkExclusionInclusion(excl1,
#                                       incl1,
#                                       excl2,
#                                       incl2) 

#       line_list[1] = e_or_i

        if hasNegativeVals(excl_raw, incl_raw, excl_lenNorm, incl_lenNorm):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            excl_raw = 0
            incl_raw = 0
            excl_lenNorm = 0
            incl_lenNorm = 0

        out_str = "%s\t%d\t%d\t%d\t%d\n" % ("\t".join(line_list), 
                                            excl_raw, incl_raw, 
                                            excl_lenNorm, incl_lenNorm)
        file2.write(out_str)

        event2counts[event_key] = (excl_raw, incl_raw, excl_lenNorm, incl_lenNorm)

    file2.close()

    return event2counts

def adjust_all_as_events(jcn_coords_str, exon_coords_str, ie_coords_str,
                         norm):

    """
    Remnant function from paired end quant
    """
    return 0
#   jcn_coords = parse_all_as_event_regions(jcn_coords_str)
#   exon_coords = parse_all_as_event_regions(exon_coords_str)
#   ie_coords = parse_all_as_event_regions(ie_coords_str)

#   overcounts = 0
#   qname_set_list = []

#   return overcounts

def allHaveSameCoord(coord_set, start_or_end):
    """
    Checks if the start or end position is the same
    """
    if start_or_end == "start":
        i = 0
    elif start_or_end == "end":
        i = 1

    this_coord_set = coord_set.copy()

    coord_pos = this_coord_set.pop()[i]

    for coord in this_coord_set:
        if coord[i] != coord_pos:
            return False

    return True

def breakInclusionRegion(alt_start_or_end, chr, ordered_pos):
    """ 
    Will divide the inclusion regions into breaks corresponding to the
    junctions.
    """
    positions = list(ordered_pos)

    excl_regions = []
    last_pos = positions.pop(0)
    if alt_start_or_end == "alt_start":
        for pos in positions:
            region = formatCoordStr(chr,
                                   last_pos,
                                   pos - 1)
            excl_regions.append(region)
            last_pos = pos
    else:
        for pos in positions:
            region = formatCoordStr(chr,
                                   last_pos + 1,
                                   pos)
            excl_regions.append(region)
            last_pos = pos

    return excl_regions
                        
            

def buildMutuallyExclusiveDict(chr, coord_start2end, coord_end2start):
    """
    Returns a dictionary of the form
    {(upstream intron start, downstream intron end):[[(exon_start, exon_end),], 
                                                     [(exon_start, exon_end),],]}
    """
    me_dict = {}
    count = 0
    for upstream_start in coord_start2end[chr]:
        count += 1

        if len(coord_start2end[chr][upstream_start]) == 1:
            continue

        for downstream_end in coord_end2start[chr]:
    
            # Reduces the search space
            if downstream_end - upstream_start + 1 > MAX_GENE_LEN:
                continue

            # Check if downstream end has more than 1 event
            if len(coord_end2start[chr][downstream_end]) == 1:
                continue

            if downstream_end < upstream_start:
                continue


            upstream_check = True
            for upstream_end in coord_start2end[chr][upstream_start]:
                if upstream_end > downstream_end:
                    upstream_check = False
                    break

            if not upstream_check:
                continue # to next downstream end

            downstream_check = True
            for downstream_start in coord_end2start[chr][downstream_end]:
                if downstream_start < upstream_start:
                    downstream_check = False
                    break

            if not downstream_check:
                continue # to next downstream end

            # these sets of introns have passed previous checks.  Now pair off
            # into exons.
            upstream_ends = list(coord_start2end[chr][upstream_start])
            downstream_starts = list(coord_end2start[chr][downstream_end])

            upstream_ends.sort()
            downstream_starts.sort()

            exon_list = []
            for i in range(len(upstream_ends)):
                for j in range(len(downstream_starts)):
                    exon_start = upstream_ends[i] + 1
                    exon_end = downstream_starts[j] - 1
                    if exon_start < exon_end:
                        if exon_end - exon_start + 1 > MAX_EXON_LEN:
                            continue
                        exon_list.append((exon_start, exon_end))

            exon_list.sort()
           
            exon_list_len = len(exon_list) 
            if exon_list_len < 2:
                continue

            if exon_list_len > MAX_EXON_CLUSTER:
#                print "MXE: Too many exons in cluster=%d" % exon_list_len
                continue

            me_exons = findNonOverlappingSets(-1, exon_list)
            me_exons = removeRedundantSets(me_exons)
        
            # Add mutually exclusive event to the dictionary
            me_list = []
            for exon_set in me_exons:
                if len(exon_set) == 1:
                    continue
                me_list.append(list(exon_set))

            # If empty set, continue
            if len(me_list) == 0:
                continue 

            me_dict[(upstream_start, downstream_end)] = me_list

    return me_dict

def checkExclusionInclusion(excl_file1_count,
                           incl_file1_count,
                           excl_file2_count,
                           incl_file2_count):
    """
    Checks the percent exclusion from both files.  
    If the percent exclusion increased from file1 to file 2, then this returns
    "E" for exclusion.
    Else, it is considered an increase in inclusion, so it will return an "I".
    
    If they are equal, it will return a "?"
    """
    # Old variable
    return "?"

#   if excl_file1_count == 0 and incl_file1_count == 0:
#       return "?"
#   if excl_file2_count == 0 and incl_file2_count == 0:
#       return "?"

#   perc_excl1 = (excl_file1_count / 
#                 (excl_file1_count + incl_file1_count))
#   perc_excl2 = (excl_file2_count / 
#                 (excl_file2_count + incl_file2_count))

#   if perc_excl2 > perc_excl1:
#       return "E"
#   elif perc_excl2 < perc_excl1:
#       return "I"

#   return "?"
    
def checkExclusionInclusion_AA_AD_AFE_ALE(alt_start_or_end,
                                          inclusion_start, inclusion_end,
                                          ordered_pos,
                                          proportions1,
                                          proportions2):
    """
    Will calculate the proportion of reads are from the inclusion isoform.  If
    the proportion is the same, then will return a ?
    
    If the proportion of the inclusion decreases in the second sample, then
    will look for junction that increases.  If the alternative splice site was
    a start position, then if junction(s) to the left increases then it was a
    shift to exclusion and if junction(s) to the right increase, then it was a
    shift to inclusion.  If both left and right junctions change, then will
    return a ?.

    Similar logic for an inclusion junction that increases or if the
    alternative splice site was an end position.
    """
    # An old variable
    return "?"

#   if sum(proportions1) == 0.0 and sum(proportions2) == 0.0:
#       return "?"

#   deltaProportions = subtract_vectors(proportions2, proportions1)

#   if alt_start_or_end == "alt_start":
#       inclusion_idx = ordered_pos.index(inclusion_start)
#   else:
#       inclusion_idx = ordered_pos.index(inclusion_end)

#   if deltaProportions[inclusion_idx] == 0.0:
#       return "?"
#  
#   # When inclusion isoform goes down, check which direction the isoform
#   # increases 
#   if deltaProportions[inclusion_idx] < 0:
#       left_increase = False
#       right_increase = False
#       
#       for i in range(0,inclusion_idx):
#           if deltaProportions[i] > 0:
#               left_increase = True
#       if inclusion_idx != (len(deltaProportions) - 1):
#           for i in range(inclusion_idx+1, len(deltaProportions)):
#               if deltaProportions[i] > 0:
#                   right_increase = True
#       if left_increase and right_increase:
#           return "?"
#       elif left_increase:
#           if alt_start_or_end == "alt_start":
#               return "E"
#           else:
#               return "I"
#       elif right_increase:
#           if alt_start_or_end == "alt_start":
#               return "I"
#           else:
#               return "E"
#       else: # No change in proportions
#           return "?"

#   # When inclusion isoform goes up, check which direction the isoform
#   # decreases
#   else: # deltaProportions[inclusion_idx] > 0.0
#       left_decrease = False
#       right_decrease = False
#       
#       for i in range(0, inclusion_idx):
#           if deltaProportions[i] < 0:
#               left_decrease = True
#       if inclusion_idx != (len(deltaProportions)-1):
#           for i in range(inclusion_idx+1,len(deltaProportions)):
#               if deltaProportions[i] < 0:
#                   right_decrease = True

#       if left_decrease and right_decrease:
#           return "?"
#       elif left_decrease:
#           if alt_start_or_end == "alt_start":
#               return "I"
#           else:
#               return "E"
#       elif right_decrease:
#           if alt_start_or_end == "alt_start":
#               return "E"
#           else:
#               return "I"
#       else: # No change in proportions
#           return "?"
    

def convertCoordStr(coord_str):
#    chr, start_str, end_str = coord_str.split()
    if ":" in coord_str:
        chr, start_end = coord_str.split(":")
        start_str, end_str = start_end.split("-")
    else:
        # Old version of coordinate string
        chr, start_str, end_str = coord_str.split("_")

    return (chr, int(start_str), int(end_str))

def coordOverlaps_wSet(this_coord, coord_set):
    for coord in coord_set:
        if coordsOverlap(this_coord[0],
                         this_coord[1],
                         coord[0],
                         coord[1]):
            return True

    return False

def copyCassetteDict(cassette_exon_dict):
    """
    Since this is a dictionary of sets, I need to make copies of the sets
    """
    new_dict = {}

    for chr in cassette_exon_dict:
        new_dict[chr] = cassette_exon_dict[chr].copy()

    return new_dict

def determineAltStartOrEnd(incl_start, incl_end, excl_jcn_list):
    """ 
    Returns whehter the start or ends are alternative
    """
    alt_start = False
    alt_end = False
    for excl_jcn in excl_jcn_list:
        chr, start, end = convertCoordStr(excl_jcn)

        if start == incl_start:
            alt_end = True
        elif end == incl_end:
            alt_start = True

    if alt_start and alt_end:
        ERROR_LOG.write("Problem in determineAltStartOrEnd\n")
        return "alt_start"

    if (not alt_start) and (not alt_end):
        ERROR_LOG.write("Problem in determineAltStartOrEnd.  No alternative.\n")
        return "alt_start"
       
    if alt_start:
        return "alt_start"

    # If not returned previously
    return "alt_end" 
            

def exonsOverlapping(chr, annotated_exons, 
                     annotated_exons_start2end,
                     annotated_exons_end2start,
                     start_or_ends, n_or_p):

    exon_list = []

    # Get all adjacent exons
    for pos in start_or_ends:
        if n_or_p == "N":
            this_exon = hasAdjExons(chr, annotated_exons_start2end, [pos], "N")
        else:
            this_exon = hasAdjExons(chr, annotated_exons_end2start, [pos], "P")
        # Only append if the adjacent exon in this group exists.
        if this_exon:
            exon_list.append(convertCoordStr(this_exon))

    # If no adjacent exons could be found, then return None.
    if exon_list == []:
        return None

    for i in range(len(exon_list)-1):
        for j in range(i,len(exon_list)):
            if not coordsOverlap(exon_list[i][1],
                                 exon_list[i][2],
                                 exon_list[j][1],
                                 exon_list[j][2]):
                return False

    return True

def filterOutExclJcns(param_all_coord_start2end, param_all_coord_end2start,
                      excl_jcns):
    """
    Returns filtered dictionaries
    """
    # Build new dictionary by filtering out exclusion junctions:
    all_coord_start2end = {}
    all_coord_end2start = {}
    

    # Build all_coord_start2end
    for chr in param_all_coord_start2end:
        for start in param_all_coord_start2end[chr]:
            for end in param_all_coord_start2end[chr][start]:
                if chr in excl_jcns:
                    if (start, end) in excl_jcns[chr]:
                        continue
                if chr in all_coord_start2end:
                    if start in all_coord_start2end[chr]:
                        all_coord_start2end[chr][start].append(end)
                    else:
                        all_coord_start2end[chr][start] = [end]
                else:
                    all_coord_start2end[chr] = {start:[end]}

    # Build all_coord_end2start
    for chr in param_all_coord_end2start:
        for end in param_all_coord_end2start[chr]:
            for start in param_all_coord_end2start[chr][end]:
                if chr in excl_jcns:
                    if (start, end) in excl_jcns[chr]:
                        continue
                if chr in all_coord_end2start:
                    if end in all_coord_end2start[chr]:
                        all_coord_end2start[chr][end].append(start)
                    else:
                        all_coord_end2start[chr][end] = [start]
                else:
                    all_coord_end2start[chr] = {end:[start]}

    return all_coord_start2end, all_coord_end2start

def findAdjacentAFE_ALE_exon(chr, afe_ale_dict, annotated_exon_dict, pos,
                             n_or_p):

    afe_ale_coord_str = hasAdjExons(chr,
                                    afe_ale_dict,
                                    [pos],
                                    n_or_p)

    if afe_ale_coord_str:
        return afe_ale_coord_str

    if not annotated_exon_dict:
        return None

    exon_coord_str = hasAdjExons(chr,
                                 annotated_exon_dict,
                                 [pos],
                                 n_or_p)

    if exon_coord_str:
        return exon_coord_str

    return None
    
def findAdjacentSharedRegion(chr, strand, annotated_exon_dict, pos, n_or_p):
    
    if chr not in annotated_exon_dict:
        return None

    if strand not in annotated_exon_dict[chr]:
        return None

    region_start_or_end = None

    for coord in annotated_exon_dict[chr][strand]:
        exon_start = coord[0]
        exon_end = coord[1]
        if exon_start <= pos <= exon_end:
            if n_or_p == "N":
                if region_start_or_end is None:
                    region_start_or_end = exon_end
                elif exon_end < region_start_or_end:
                    region_start_or_end = exon_end
            else: # n_or_p == "P"
                if region_start_or_end is None:
                    region_start_or_end = exon_start
                elif exon_start > region_start_or_end:
                    region_start_or_end = exon_start

    if region_start_or_end is None:
        return None

    if n_or_p == "N":
        return formatCoordStr(chr, pos, region_start_or_end)

    return formatCoordStr(chr, region_start_or_end, pos)


def findExclSet(prox_excl_jcn, excl_jcns_full):
    """ 
    Returns the comma delimited list of exclusion junctions that the proximal
    junction is grouped with. Used to udpate AFE/ALE exon counts
    """
    for excl_set in excl_jcns_full:
        excl_jcns = excl_set.split(",")
    
        if prox_excl_jcn in excl_jcns:
            return excl_set


    # If not returned, there is some error
    print "Error with AFE/ALE exclusion exons: %s" % excl_jcns_full
    sys.exit(1)

    return None

def findNonOverlappingSets(current_pos, coord_list):
    """
    Recursive function.  Given the input coord_list, the function will return
    all possible non_overlapping sets of coordinates.
    Coord list will be a list of (start, end) tuples.
    """
    start_idx = 0
    end_idx = 1
    new_set = []

    for coord in coord_list:
        if coord[start_idx] > current_pos:
            new_set.append(set([coord]))
    
            additional_sets = findNonOverlappingSets(coord[end_idx],
                                                     coord_list)

            for set1 in additional_sets:
                set1.add(coord)
                new_set.append(set1)

    return new_set

def fixIRExclusion_count(ir_file_name, all_jcn_count_dict,
                         norm1, norm2, jcn_seq_len):

    file = open(ir_file_name)
    lines = file.readlines()
    file.close()

    file2 = open(ir_file_name, "w")

    for line in lines:
        line = line.rstrip("\n")

        line_list = line.split("\t")

        excl_jcn_str = line_list[5]

#        jcn_ct1 = all_jcn_count_dict[excl_jcn_str][0]
        jcn_ct_raw = all_jcn_count_dict[excl_jcn_str][1]

        if norm2:
#            jcn_ct1 = int(round(jcn_ct1/norm1))
            jcn_ct_raw = int(round(jcn_ct_raw/norm2))

#        jcn_ct1 = normalizeByLen(jcn_ct1, jcn_seq_len)
        jcn_ct_lenNorm = normalizeByLen(jcn_ct_raw, jcn_seq_len)

        if hasNegativeVals(jcn_ct_raw, jcn_ct_lenNorm, 0, 0):
            ERROR_LOG.write("Negative Vals: %s\n" % line)
            jcn_ct_raw = 0
            jcn_ct_lenNorm = 0
    
        line_list[11] = repr(jcn_ct_raw)

        line_list[13] = repr(jcn_ct_raw)

        outline = "\t".join(line_list)
        file2.write(outline + "\n")

    file2.close()

def formatChr(chr):
    """
    The chromosome string must begin with "chr"
    """
    if not chr.startswith("chr"):
        return "chr" + chr
    
    return chr

def formatCoordStr(chr, start, end):
    try:
        return "%s:%d-%d" % (chr, start, end)
    except:
        pass
        
    return "%s:%d-%d" % (chr,int(start), int(end))

    

def formatLine(line):
    line = line.strip()
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getAllEventStr(label, type, gene_name, chr, strand, 
                   excl_jcns, incl_jcns,
                   excl_exons, incl_exons,
                   incl_i_e_jcns,
                   const_exons,
                   excl_jcn_counts_raw, incl_jcn_counts_raw,
                   sum_excl_jcns_raw, sum_incl_jcns_raw,
                   excl_exon_counts_raw, incl_exon_counts_raw,
                   sum_excl_exon_raw, sum_incl_exon_raw,
                   ie_jcn_cts_raw, 
                   sum_ie_jcn_cts_raw, 
                   const_exon_cts_raw,
                   sum_const_exon_cts_raw):

    if sum_excl_jcns_raw is None:
        sum_excl_jcn_raw = ""
    if sum_incl_jcns_raw is None:
        sum_incl_jcns_raw = ""

    if sum_excl_exon_raw is None:
        sum_excl_exon_raw = ""
    if sum_incl_exon_raw is None:
        sum_incl_exon_raw = ""

    if sum_ie_jcn_cts_raw is None:
        sum_ie_jcn_cts_raw = "" 
    
    if sum_const_exon_cts_raw is None:
        sum_const_exon_cts_raw = "" 
                   
    out_str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (label,
                                                                          type,
                                                                          gene_name,
                                                                          chr,
                                                                          strand,
                                                                          excl_jcns, incl_jcns,
                                                                          excl_exons, incl_exons,
                                                                          incl_i_e_jcns,
                                                                          const_exons,
                                                                          excl_jcn_counts_raw, incl_jcn_counts_raw,
                                                                          sum_excl_jcns_raw,
                                                                          sum_incl_jcns_raw)
    out_str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (excl_exon_counts_raw, 
                   incl_exon_counts_raw,
                   sum_excl_exon_raw, 
                   sum_incl_exon_raw,
                   ie_jcn_cts_raw, 
                   sum_ie_jcn_cts_raw, 
                   const_exon_cts_raw, 
                   sum_const_exon_cts_raw)

    return out_str


def getFirstLastExons(db, txt_db, this_chr):
    """
    Returns coordinates for all annotated alternative initiation and
    termination exons

    Returns first exon dict, last exon dict in the form
    # {chr:(start, end)} 

    (alt_first_exons, 
     alt_last_exons,
     alt_first_exons_start2end,
     alt_first_exons_end2start,
     alt_last_exons_start2end,
     alt_last_exons_end2start)
    """
    afe_dict = {}
    ale_dict = {}

    afe_dict_start2end = {}
    afe_dict_end2start = {}

    ale_dict_start2end = {}
    ale_dict_end2start = {}

    # {chr:coordSearchTree}
    afe_search_tree = {}
    ale_search_tree = {}

    # {gene:set([(exon_start, exon_end)])}
    gene2most_upstrm_exons = {}
    gene2most_dwnstrm_exons = {}

    # {txt:set([(exon_start, exon_end)])
    transcript2exons = {}

    # {gene:(chr, strand)
    gene2chr_strand = {}

    # {gene: set([txt_id])}
    gene2transcripts = {}

    # First get all transcripts:
    txt_select = """SELECT DISTINCT transcript_id, gene_name, chr,
                                   strand, start, end
                          FROM exon"""
   
#    txt_select = """SELECT transcript_id, gene_name, chr,
#                           strand, start, end
#                    FROM exon"""

    if this_chr:
        txt_select += " WHERE chr = \'%s\' OR chr = \'%s\'" % (this_chr,
                                                                 this_chr.lstrip("chr"))

    txt_records = db.getDBRecords_Dict(txt_select, txt_db) 

    all_chr_set = set([])

    # Populate dictionaries        
    for exon_row in txt_records:
        txt_id = exon_row["transcript_id"]
        gene_name = exon_row["gene_name"]
        chr = formatChr(exon_row["chr"])
        strand = exon_row["strand"]

        exon_start = int(exon_row["start"])
        exon_end = int(exon_row["end"])

        if strand != "+" and strand != "-":
            print "Can't determine strand for gene: %s" % gene_name
            continue

        all_chr_set.add(chr)

        if gene_name in gene2chr_strand:
            if (chr, strand) != gene2chr_strand[gene_name]:
                print "Gene information does not agree in db: %s" % gene_name
                continue
        else:
            gene2chr_strand[gene_name] = (chr, strand)
    
        updateDictOfSets(gene2transcripts, gene_name, txt_id)
        updateDictOfSets(transcript2exons, txt_id, (exon_start, exon_end))
    
    # For all transcripts of the same gene, get the most upstream and downstrm
    # exons
    for gene in gene2transcripts:
        this_chr = gene2chr_strand[gene][0]
        for txt in gene2transcripts[gene]:
            exon_list = list(transcript2exons[txt])

            # Skip single exon genes
            if len(exon_list) == 1:
                continue

            exon_list.sort()

#           updateDictOfSets(gene2most_upstrm_exons, gene, exon_list[0])
#           updateDictOfSets(gene2most_dwnstrm_exons, gene, exon_list[-1])

            
            if gene2chr_strand[gene][1] == "+":
                updateDictOfLists(afe_dict, this_chr, exon_list[0])
                updateDictOfLists(ale_dict, this_chr, exon_list[-1])

                updateDictExons(afe_dict_start2end, this_chr, 
                                exon_list[0][0], exon_list[0][1])
                updateDictExons(afe_dict_end2start, this_chr, 
                                exon_list[0][1], exon_list[0][0])

                updateDictExons(ale_dict_start2end, this_chr, 
                                exon_list[-1][0], exon_list[-1][1])
                updateDictExons(ale_dict_end2start, this_chr, 
                                exon_list[-1][1], exon_list[-1][0])


            else: 
                updateDictOfLists(afe_dict, this_chr, exon_list[-1])
                updateDictOfLists(ale_dict, this_chr, exon_list[0])

                updateDictExons(afe_dict_start2end, this_chr, 
                                exon_list[-1][0], exon_list[-1][1])
                updateDictExons(afe_dict_end2start, this_chr, 
                                exon_list[-1][1], exon_list[-1][0])

                updateDictExons(ale_dict_start2end, this_chr, 
                                exon_list[0][0], exon_list[0][1])
                updateDictExons(ale_dict_end2start, this_chr, 
                                exon_list[0][1], exon_list[0][0])


    # For every gene, look for a unique end or start position for upstrm and
    # dwnstrm exons
    # 120228: No longer enforcing that the first and last exons be alternative
    # since this is queried in printAlternativeDonorsAcceptors.
#   for gene in gene2most_upstrm_exons:
#       exon_list = list(gene2most_upstrm_exons[gene])
#       if len(exon_list) == 1:
#           continue
#       this_exon_end = exon_list[0][1]
#       hasDifferentEnds = False
#       for exon in exon_list:
#           if this_exon_end != exon[1]:
#               hasDifferentEnds = True
#       if hasDifferentEnds:
#           # If strand is + then it is an AFE
#           if gene2chr_strand[gene][1] == "+":
#               for exon in exon_list:
#                   updateDictOfLists(afe_dict, gene2chr_strand[gene][0], exon)
#           else: 
#               for exon in exon_list:
#                   updateDictOfLists(ale_dict, gene2chr_strand[gene][0], exon)

#   for gene in gene2most_dwnstrm_exons:
#       exon_list = list(gene2most_dwnstrm_exons[gene])
#       if len(exon_list) == 1:
#           continue
#       this_exon_start = exon_list[0][0]
#       hasDifferentStarts = False
#       for exon in exon_list:
#           if this_exon_start != exon[0]:
#               hasDifferentStarts = True
#       if hasDifferentStarts:
#           # If strand is + then it is an ALE
#           if gene2chr_strand[gene][1] == "+":
#               for exon in exon_list:
#                   updateDictOfLists(ale_dict, gene2chr_strand[gene][0], exon)
#           else: 
#               for exon in exon_list:
#                   updateDictOfLists(afe_dict, gene2chr_strand[gene][0], exon)

    # Add chromsomes that did not have AFE or ALE to avoid key errors
    for chr in all_chr_set:
        if chr not in afe_dict:
            afe_dict[chr] = []

        if chr not in ale_dict:
            ale_dict[chr] = []

        if chr not in afe_dict_start2end:
            afe_dict_start2end[chr] = {}

        if chr not in afe_dict_end2start:
            afe_dict_end2start[chr] = {}

        if chr not in ale_dict_start2end:
            ale_dict_start2end[chr] = {}

        if chr not in ale_dict_end2start:
            ale_dict_end2start[chr] = {}

    for chr in afe_dict:
        afe_search_tree[chr] = getSearchTree(afe_dict[chr])
    for chr in ale_dict:
        ale_search_tree[chr] = getSearchTree(ale_dict[chr])

    return (afe_dict, ale_dict, 
            afe_dict_start2end, afe_dict_end2start,
            ale_dict_start2end, ale_dict_end2start,
            afe_search_tree, ale_search_tree)

def updateDictExons(exon_dict, chr, pos1, pos2):
    try:
        updateDictOfSets(exon_dict[chr], pos1, pos2)
    except:
        exon_dict[chr] = {pos1:set([pos2])}

def getAnnotatedExonCoords(db, txt_db, this_chr):
    """
    Returns a dictionary of {chr:set([(exon_start, exon_end)])}
    """
    exon_select = """SELECT DISTINCT transcript_id, chr, strand, start, end
                       FROM exon"""

    if this_chr:
        exon_select += " WHERE chr = \'%s\' OR chr = \'%s\'" % (this_chr,
                                                                this_chr.lstrip("chr"))
    exon_records = db.getDBRecords_Dict(exon_select, txt_db)

    exon_dict = {}
    exon_internal_dict = {}
    exon_dict_no_strand = {}
    exon_dict_by_strand = {}
    
    exon_search_tree_coords = {}
    

    transcript2exons = {}

    for exon_row in exon_records:
        txt_id = exon_row["transcript_id"]
        chr = formatChr(exon_row["chr"])
        start = int(exon_row["start"])
        end = int(exon_row["end"])
        strand = exon_row["strand"]

        if chr in exon_dict:
            exon_dict[chr].add((start, end, strand))
            exon_dict_no_strand[chr].add((start, end))
            try:
                exon_dict_by_strand[chr][strand].add((start,end))
            except:
                exon_dict_by_strand[chr] = {strand: set([(start,end)])}
        else:
            exon_dict[chr] = set([(start, end, strand)])
            exon_dict_no_strand[chr] = set([(start, end)])
            exon_dict_by_strand[chr] = {strand:set([(start,end)])}

        updateDictOfSets(transcript2exons, txt_id, (chr, start, end, strand))

        # Add to search tree
        if chr in exon_search_tree_coords:
            if strand in exon_search_tree_coords[chr]:
                exon_search_tree_coords[chr][strand].add((start,end))
            else:
                exon_search_tree_coords[chr][strand] = set([(start,end)])
        else:
            exon_search_tree_coords[chr] = {strand:
                                            set([(start,end)])}                                      

    for txt_id in transcript2exons:
        exon_list = []

        chr = None

        for exon_coord in transcript2exons[txt_id]:
            chr = exon_coord[0]
            exon_list.append((exon_coord[1],
                              exon_coord[2]))
        exon_list.sort()

        # Only add internal exons
        for i in range(1, len(exon_list)-1):
            updateDictOfSets(exon_internal_dict, chr, exon_list[i])

    # Make search tree
    exon_search_tree = {}
    for chr in exon_search_tree_coords:
        exon_search_tree[chr] = {}
        for strand in exon_search_tree_coords[chr]:
            exon_search_tree[chr][strand] = getSearchTree(list(exon_search_tree_coords[chr][strand]))

        # Add empty search trees for remaining strands
        for strand in ["+", "-", "."]:
            if strand not in exon_search_tree_coords[chr]:
                exon_search_tree[chr][strand] = getSearchTree([])

    return exon_dict, exon_internal_dict, exon_dict_no_strand, exon_dict_by_strand, exon_search_tree

    # {chr: (start,end): [gene_names]}
def getAnnotatedGenes(db, txt_db, this_chr):
    intron_select = """SELECT gene_name, chr, strand,
                            start, end
                     FROM intron"""

    if this_chr:
        intron_select += " WHERE chr = \'%s\' OR chr = \'%s\'" % (this_chr,
                                                                this_chr.lstrip("chr"))

    intron_records = list(db.getDBRecords_Dict(intron_select, txt_db))

    exon_select = """SELECT gene_name, chr, strand, start, end FROM exon"""
    if this_chr:
        exon_select += " WHERE chr = \'%s\' OR chr = \'%s\'" % (this_chr,
                                                                  this_chr.lstrip("chr"))

    exon_records = list(db.getDBRecords_Dict(exon_select, txt_db))

    all_records = intron_records + exon_records

    gene_dict = {}
    gene_dict_by_strand = {}
    for row in all_records:

#        gid = int(row["gene_name"])
        name = row["gene_name"]

        strand = row["strand"]
        start = int(row["start"])
        end = int(row["end"])

        chr = formatChr(row["chr"])

        try:
            if (start,end) in gene_dict[chr]:
                updateDictOfLists(gene_dict[chr], (start,end), name)
            else:
                gene_dict[chr][(start,end)] = [name]
#            updateDictOfLists(gene_dict[chr], chr, (gid, name, strand, start, end)) 
        except:
            gene_dict[chr] = {(start,end):[name]}

        try:
            if (start,end) in gene_dict_by_strand[chr][strand]:
                updateDictOfLists(gene_dict_by_strand[chr][strand], (start,end), name)
            else:
                gene_dict_by_strand[chr][strand][(start,end)] = [name]
        except:
            if chr in gene_dict_by_strand:
                if strand in gene_dict_by_strand:
                    updateDictOfLists(gene_dict_by_strand[chr][strand], (start,end), name) 
                else:
                    gene_dict_by_strand[chr][strand] = {(start,end):[name]}
            else:
                gene_dict_by_strand[chr] = {strand:{(start,end):[name]}}

#       if chr in gene_dict_by_strand:
#           updateDictOfLists(gene_dict_by_strand[chr], strand, (gid, name, strand, start, end))
#       else:
#           gene_dict_by_strand[chr] = {strand:
#                                       [(gid, name, strand, start,end)]}

    # If for some reason there are no genes on a chromosome, can check for
    # chromosomes that have exon to avoid key errors on chromosomes
#   exon_select = "SELECT DISTINCT chr FROM exon"

#   exon_records = db.getDBRecords_Dict(exon_select, txt_db)

#   for row in exon_records:
#       chr = formatChr(row["chr"])

#       if chr not in gene_dict:
#           # Add an empty list.
#           gene_dict[chr] = []

#       # Update gene_dict_by_strand as well
#       if chr not in gene_dict_by_strand:
#           gene_dict_by_strand[chr] = {"+": [],
#                                       "-": []}
#       else:
#           if "+" not in gene_dict_by_strand[chr]:
#               gene_dict_by_strand[chr]["+"] = []
#           if "-" not in gene_dict_by_strand[chr]:
#               gene_dict_by_strand[chr]["-"] = [] 

    return gene_dict, gene_dict_by_strand

def getAnnotatedIntronCoords(db, txt_db, this_chr):
    """
    Returns a dictionary of {chr:set([(intron_start, intron_end)])}
    """
    intron_select = """SELECT DISTINCT chr, strand, start, end
                       FROM intron"""

    if this_chr:
        intron_select += " WHERE chr = \'%s\' OR chr = \'%s\'" % (this_chr,
                                                                  this_chr.lstrip("chr"))

    intron_records = db.getDBRecords_Dict(intron_select, txt_db)

    intron_dict = {}

    for intron_row in intron_records:
        chr = formatChr(intron_row["chr"])
        start = int(intron_row["start"])
        end = int(intron_row["end"])
        strand = intron_row["strand"]

        if chr in intron_dict:
            intron_dict[chr].add((start, end, strand))
        else:
            intron_dict[chr] = set([(start, end, strand)])

    # Add chromsomes that do not contain introns to the dictionary
    exon_select = "SELECT DISTINCT chr FROM exon"

    exon_records = db.getDBRecords_Dict(exon_select, txt_db)

    for exon_row in exon_records:
        chr = formatChr(exon_row["chr"])

        if chr not in intron_dict:
            # add empty set since there are no introns on this chromosome
            intron_dict[chr] = set([])

    return intron_dict
    

def getColSums(cols, line_elems):
    sum = 0

    for col in cols:
        if line_elems[col] != "":
            ct = int(line_elems[col])
            sum += ct

    return sum

def getCoordCounts4all_as_events(all_as_exon_str, 
                                 mapped_file1_counts, mapped_file2_counts):
    """ 
    Parses the exon string, finds the corresponding counts, then outputs counts
    in proper format
    Outputs
    exon_ct_str_samp1, exon_ct_str_samp2,
    sum_exon_ct_samp1, sum_exon_ct_samp2
    """
    if "," in all_as_exon_str:
        print "Need to account for , in exon_str list"
        sys.exit(1)

    if all_as_exon_str == "" or all_as_exon_str == "None":
        return "", "", None, None

    exon_list = all_as_exon_str.split(";")

    sum_samp1 = 0
    samp1_cts = []
    sum_samp2 = 0
    samp2_cts = []

    for exon_str in exon_list: 
        try:
            samp1_ct = mapped_file1_counts[exon_str]
            samp2_ct = mapped_file2_counts[exon_str]

        except:
            ERROR_LOG.write("Could not find counts for %s\n" % exon_str)
            samp1_ct = 0
            samp2_ct = 0

        sum_samp1 += samp1_ct
        samp1_cts.append(samp1_ct)

        sum_samp2 += samp2_ct
        samp2_cts.append(samp2_ct)

    return (";".join(map(repr, samp1_cts)), ";".join(map(repr, samp2_cts)),
            sum_samp1, sum_samp2)

def getExonDistanceDifference(exon_str_list):
    """
    Takes a list of exon strings and returns the maximum distance between any
    pair of the exons. 
    """
    max_exon_len = 0
    min_exon_len = INFINITY

    shortest_exon_str = ""

    for exon_str in exon_str_list:
        chr, start, end = convertCoordStr(exon_str)
#       chr, start_str, end_str = exon_str.split()

#       start = int(start_str)
#       end = int(end_str)
        
        exon_len = end - start + 1

        if exon_len > max_exon_len:
            max_exon_len = exon_len
        if exon_len < min_exon_len:
            min_exon_len = exon_len
            shortest_exon_str = exon_str

    return max_exon_len - min_exon_len, shortest_exon_str

    
def getInclusionPortion(inclusion_exon, other_exon_list):
    """
    Input is an inclusion exon, and a list of other related exons.
    Returns inclusion_start, inclusion_end
    """
    chr, inclusion_start, inclusion_end = convertCoordStr(inclusion_exon)
#   chr, inclusion_start_str, inclusion_end_str = inclusion_exon.split()

#   inclusion_start = int(inclusion_start_str)
#   inclusion_end = int(inclusion_end_str)

    for exon_str in other_exon_list:
        chr, other_start_str, other_end_str = convertCoordStr(exon_str)
#       chr, other_start_str, other_end_str = exon_str.split)
#  
#       other_start = int(other_start_str)
#       other_end = int(other_end_str) 
   
        if inclusion_start < other_end < inclusion_end: 
            inclusion_start = other_end + 1
            continue
        elif inclusion_start < other_start < inclusion_end:
            inclusion_end = other_start - 1
    
    return inclusion_start, inclusion_end
    
    
def getInternalIntrons(intron_search_tree, start_pos, end_pos):
    """
    Will return a list of all intron coordinates that are contained within the
    start and end positions
    """
    intron_coords = []
    
    this_coord = (start_pos, end_pos)

    findInternalCoords(intron_search_tree, intron_coords, this_coord)    
#   for intron_start in coord_start2end[chr]:
#       if intron_start > start_pos:
#           for intron_end in coord_start2end[chr][intron_start]:
#               if intron_end < end_pos:
#                   intron_coords.append((intron_start, intron_end))

    if intron_coords == []:
        return None

    return intron_coords

def getJcnStrandInfo(jcn2strand, genome_file):
    """
    Organizes junctions into a new dictionary that is keyed by chromosome, then
    opens the chromosome files to get the junction sequences.
    """
    print "Inferring Junction Strand Info"
    chr2jcns = {}
    for jcn in jcn2strand:

        # Only check for strand information if it does not exist
        if jcn2strand[jcn] == "+" or jcn2strand[jcn] == "-":
            continue

        (chr, start, end) = convertCoordStr(jcn)     
   
        if chr in chr2jcns:
            chr2jcns[chr].append((chr, start, end))
        else:
            chr2jcns[chr] = [(chr, start, end)] 

    # Now opening genome sequence
    try:
        genome_fileh = open(genome_file)
    except:
        print "Could not find the genome sequence."
        sys.exit(1)

    for record in SeqIO.parse(genome_fileh, "fasta"):
        file_chr = formatChr(record.id)
        chr_seq = str(record.seq)

        print "Inferring junctions on chromosome: %s" % file_chr

        if file_chr in chr2jcns:
            for (chr, start, end) in chr2jcns[file_chr]:

                jcn_str = formatCoordStr(chr, start, end)
                intron_seq = chr_seq[start-1:end]

                if intron_seq.startswith("GT") and intron_seq.endswith("AG"):
                    jcn2strand[jcn_str] = "+"
                elif intron_seq.startswith("CT") and intron_seq.endswith("AC"):
                    jcn2strand[jcn_str] = "-"
                # Other common splice site sequence
                elif intron_seq.startswith("GC") and intron_seq.endswith("AG"):
                    jcn2strand[jcn_str] = "+"
                elif intron_seq.startswith("CT") and intron_seq.endswith("GC"):
                    jcn2strand[jcn_str] = "-"
                # minor spliceosome
                elif intron_seq.startswith("AT") and intron_seq.endswith("AC"):
                    jcn2strand[jcn_str] = "+"
                elif intron_seq.startswith("GT") and intron_seq.endswith("AT"):
                    jcn2strand[jcn_str] = "-"
                # Priority to 5' splice site since there is more information
                # there
                elif intron_seq.startswith("GT"):
                    jcn2strand[jcn_str] = "+"
                elif intron_seq.endswith("AC"):
                    jcn2strand[jcn_str] = "-"
                elif intron_seq.endswith("AG"):
                    jcn2strand[jcn_str] = "+"
                elif intron_seq.startswith("CT"):
                    jcn2strand[jcn_str] = "-"
                else:
                    print "Cannot find strand for %s" % jcn_str
                    
def getJunctionProportion(chr, isoform_jcn,
                          all_jcn_count_dict,
                          anchor_type,
                          jcn_end_dict,
                          anchor_pos):
    """ 
    Pseudo counts are added to prevent proportions that are 0
    """
    isoform_jcn_ct1 = all_jcn_count_dict[isoform_jcn][0] + 1
    isoform_jcn_ct2 = all_jcn_count_dict[isoform_jcn][1] + 1

    total_ct1 = 0
    total_ct2 = 0

    for other_pos in jcn_end_dict[chr][anchor_pos]:

        if anchor_type == "start":
            this_jcn = formatCoordStr(chr, anchor_pos, other_pos)
        else: # anchor_type == end
            this_jcn = formatCoordStr(chr, other_pos, anchor_pos)

        total_ct1 += all_jcn_count_dict[this_jcn][0] + 1 
        total_ct2 += all_jcn_count_dict[this_jcn][1] + 1 

    
    return (float(isoform_jcn_ct1)/total_ct1, float(isoform_jcn_ct2)/total_ct2)


def getJunctionProportion_AFE_ALE(chr, isoform_jcns,
                                  all_jcn_count_dict,
                                  anchor_type,
                                  jcn_end_dict):
    """ 
    Pseudo counts are added to prevent proportions that are 0
    """
    isoform_jcns_list = isoform_jcns.split(",")

    isoform_jcn_ct1 = 0
    isoform_jcn_ct2 = 0

    total_ct1 = 0
    total_ct2 = 0

    anchor_positions = set([])

    for isoform_jcn in isoform_jcns_list:
        isoform_jcn_ct1 += all_jcn_count_dict[isoform_jcn][0] + 1
        isoform_jcn_ct2 += all_jcn_count_dict[isoform_jcn][1] + 1

        this_chr, this_start, this_end = convertCoordStr(isoform_jcn)
    
        if anchor_type == "start":
            anchor_positions.add(this_start)         
        else: # end
            anchor_positions.add(this_end)

    for anchor_pos in anchor_positions:
        for other_pos in jcn_end_dict[chr][anchor_pos]:

            if anchor_type == "start":
                this_jcn = formatCoordStr(chr, anchor_pos, other_pos)
            else: # anchor_type == end
                this_jcn = formatCoordStr(chr, other_pos, anchor_pos)

            total_ct1 += all_jcn_count_dict[this_jcn][0] + 1 
            total_ct2 += all_jcn_count_dict[this_jcn][1] + 1 

    
    return (float(isoform_jcn_ct1)/total_ct1, float(isoform_jcn_ct2)/total_ct2)

def getLongestJunction(junctions):
    """
    From a list of junctions get the longest
    """
    longest_jcn = None
    longest_len = -1

    for jcn in junctions:
        chr, start, end = convertCoordStr(jcn)
        if (end - start + 1) > longest_len:
            longest_jcn = jcn
            longest_len = end - start + 1

    return longest_jcn

def getMXE_MCProportion(exon_coord, upstrm_jcn_start, dwnstrm_jcn_end,
                        all_jcn_count_dict, 
                        all_coord_start2end,
                        all_coord_end2start):
    """ 
    Returns the proportion of junction reads that support the inclusion of this
    mutually exclusive exon or multi-cassette exon in both samples.
    """
    chr, exon_start, exon_end = convertCoordStr(exon_coord)
    upstrm_jcn = formatCoordStr(chr, upstrm_jcn_start, exon_start - 1)
    dwnstrm_jcn = formatCoordStr(chr, exon_end + 1, dwnstrm_jcn_end)

    (upstrm_prop1,
     upstrm_prop2) = getJunctionProportion(chr,
                                           upstrm_jcn,
                                           all_jcn_count_dict,
                                           "end",
                                           all_coord_end2start,
                                           exon_start - 1)
       

    (dwnstrm_prop1,
     dwnstrm_prop2) = getJunctionProportion(chr,
                                            dwnstrm_jcn,
                                            all_jcn_count_dict,
                                            "start",
                                            all_coord_start2end,
                                            exon_end + 1)
 

    # Take the minimum proportion from both sides
    return (min(upstrm_prop1, dwnstrm_prop1), 
            min(upstrm_prop2, dwnstrm_prop2))
        


def getPairedCoordAdjust(qname_set_list, norm):
    """
    Input is a list of sets.  The list has an arbitrary length.  This function
    count the 
    """
    # Any qname with more than a count of 1 will flag as overcounting
    qname2count = {}

    for setlist in qname_set_list:
        for qname in setlist:
            if qname in qname2count:
                qname2count[qname] += 1
            else:
                qname2count[qname] = 1

    multiple_count = 0
    for qname in qname2count:
        if qname2count[qname] > 1:
            overcount = qname2count[qname] - 1 
            multiple_count += overcount

    if norm:
        multiple_count = int(round(multiple_count/norm))

    return multiple_count

def getPairedRegionCoordAdjust(paired_coord2qname2count, region_coord, norm):
    """ 
    Returns the number of duplicated counts.
    """
    multiple_count = 0
    if region_coord in paired_coord2qname2count:
        for qname in paired_coord2qname2count[region_coord]: 
            q_count = paired_coord2qname2count[region_coord][qname]

            if q_count > 1:
                overcount = q_count - 1
                multiple_count += overcount

    if norm:
        multiple_count = int(round(multple_count/norm))

    return multiple_count

def getPossibleLeftEnds(possible_excl, excl_end):
    """
    Returns list of any values less than the index of the excl_end
    """
    return list(possible_excl[:possible_excl.index(excl_end)])

def getSSOrderAndProportions(alt_start_or_end, inclusion_start, inclusion_end,
                             exclusion_str_list,
                             exclusion_cts1_list,
                             inclusion_cts1,
                             exclusion_cts2_list,
                             inclusion_cts2):
    """ 
    Returns three lists:
    ordered_pos - the splice site positions in order
    proportions1 - the proportion of each isoform in sample 1
    proportions2 - the proportion of each isoform in sample 2
    ordered_cts1_list
    ordered_cts2_list
    """
    unordered_pos = []
    unordered_prop1 = []
    unordered_prop2 = []

    unordered_inclusion_pos = None

    if alt_start_or_end == "alt_start":
        unordered_pos.append(inclusion_start)
        unordered_inclusion_pos = inclusion_start
        for exclusion_str in exclusion_str_list:
            chr, start, end = convertCoordStr(exclusion_str)
            unordered_pos.append(start)
    else:
        unordered_pos.append(inclusion_end)
        unordered_inclusion_pos = inclusion_end
        for exclusion_str in exclusion_str_list:
            chr, start, end = convertCoordStr(exclusion_str)
            unordered_pos.append(end)

    unordered_prop1.append(inclusion_cts1)
    for excl_cts in exclusion_cts1_list:
        unordered_prop1.append(excl_cts)
    
    unordered_prop2.append(inclusion_cts2)
    for excl_cts in exclusion_cts2_list:
        unordered_prop2.append(excl_cts)

    sum_raw = sum(unordered_prop1)
    sum_lenNorm = sum(unordered_prop2)

    ordered_pos = list(unordered_pos)
    ordered_pos.sort()

    ordered_prop1 = []
    ordered_prop2 = []
    ordered_cts1 = []
    ordered_cts2 = []

    num_isoforms = len(ordered_pos)
    
    for pos in ordered_pos:
        unordered_idx = unordered_pos.index(pos)

        if sum_raw == 0:
            ordered_prop1.append(0.0)
            ordered_cts1.append(0)
        else:
            # Add pseudo count for proportions only, not real counts
            ordered_prop1.append(float((unordered_prop1[unordered_idx]+1)/(sum_raw+num_isoforms)))
            ordered_cts1.append(unordered_prop1[unordered_idx])
            
        
        if sum_lenNorm == 0:
            ordered_prop2.append(0.0)
            ordered_cts2.append(0)
        else:
            # Add pseudo count for proportions only, not real counts
            ordered_prop2.append(float((unordered_prop2[unordered_idx]+1)/(sum_lenNorm+num_isoforms)))
            ordered_cts2.append(unordered_prop2[unordered_idx])
               
    return (ordered_pos, ordered_prop1, ordered_prop2, ordered_cts1, ordered_cts2)

def hasAdjExons(chr, exon_dict, start_or_ends, n_or_p):
    """ 
    Looks for adjacent alternative exons next or previous to the intron start and
    ends that are given.
    """
    if chr not in exon_dict:
        return None

    # If there are multiple adjacent exons, pick the longest
    longest_upstrm_start = INFINITY
    longest_dwnstrm_end = -1

    for pos in start_or_ends:
#        for (exon_start, exon_end) in exon_dict[chr]:
        if n_or_p == "N":
            exon_start = pos + 1
            if exon_start in exon_dict[chr]:
                for exon_end in exon_dict[chr][exon_start]:
                    if exon_end > longest_dwnstrm_end:
                        longest_upstrm_start = exon_start
                        longest_dwnstrm_end = exon_end
        else: # n_or_p == "P"
            exon_end = pos - 1
            if exon_end in exon_dict[chr]:
                for exon_start in exon_dict[chr][exon_end]:
                    if exon_start < longest_upstrm_start:
                        longest_upstrm_start = exon_start
                        longest_dwnstrm_end = exon_end

    if longest_dwnstrm_end != -1:
        return formatCoordStr(chr, longest_upstrm_start, longest_dwnstrm_end)

    # Was not returned previously
    return None 
    


def hasInternalIntron(coord_start2end, chr, exon_start, exon_end):

    for intron_start in coord_start2end[chr]:
        if intron_start >= exon_start and intron_start <= exon_end:
            for intron_end in coord_start2end[chr][intron_start]:
                if intron_end >= exon_start and intron_end <= exon_end:
                    return True

    return False

def hasInternalTerminationInitationExon(alt_exon_dict, chr, exon_start, exon_end):
    """
    I noticed that some novel events ended up being termination or initiation
    exons.
    """
    if chr not in alt_exon_dict:
        return False

    for start, end in alt_exon_dict[chr]:
        if exon_start == start:
            return True
        if exon_end == end:
            return True

    return False

def hasNegativeVals(excl1, incl1, excl2, incl2):
    if excl1 < 0:
        return True
    if incl1 < 0:
        return True
    if excl2 < 0:
        return True
    if incl2 < 0 :
        return True

    return False

def hasOtherOverlappingGene(annotated_genes, alt_polyA_dict, chr, sgid, anchor):
    """
    Will return true if there is an overlapping gene but on a different strand
    """
    # Find largest exon in group
    left_most_coord = INFINITY
    right_most_coord = 0
   
    dict_key = (chr, sgid, anchor) 
   
    for exon_str in alt_polyA_dict[dict_key]:
        exon_chr, start_str, end_str = convertCoordStr(exon_str)
#       exon_chr, start_str, end_str = exon_str.split() 

#       start = int(start_str)
#       end = int(end_str)

        if start < left_most_coord:
            left_most_coord = start
        if end > right_most_coord:
            right_most_coord = end

    # annotated_genes -  {chr: [(sgid, name, strand, start, end)]}
    for gene_info in annotated_genes[chr]:
        if gene_info[0] == sgid:
            continue

        if coordsOverlap(left_most_coord, right_most_coord,
                         gene_info[3], gene_info[4]):
            return True


    return False   
 
def inferExclusionJunctions(upstrm_jcns, dwnstrm_jcns):
    excl_jcns = []

    for upstrm_jcn in upstrm_jcns:
        chr, upstr_start, upstr_end = convertCoordStr(upstrm_jcn)
        for dwnstr_jcn in dwnstrm_jcns:
            cnr, dwnstr_start, dwnstr_end = convertCoordStr(dwnstr_jcn)

            excl_jcn = formatCoordStr(chr, upstr_start, dwnstr_end)
            excl_jcns.append(excl_jcn)
   
    return excl_jcns 

def inferGeneName(annotated_genes_by_strand, chr, strand, start_and_ends):
    """
    Looks for the gene(s) that the start and end coordinate is contained in
    """
    # annotated_genes -  
    # {chr: (start,end): [gene_names]}
    # {chr: strand: (start,end): [gene_names]}
    gene_set = set([])
    if chr not in annotated_genes_by_strand:
        return "None"

    if strand not in annotated_genes_by_strand[chr]:
        return "None"

    for elem1 in start_and_ends:
        elem2 = []
        if "," in elem1:
            elem2 = elem1.split(",")
        else:
            elem2 = [elem1]
        for elem in elem2:
            if ":" in elem:
                chr, start, end = convertCoordStr(elem)
                elem = (start, end)

            if elem in annotated_genes_by_strand[chr][strand]:
                gene_set.update(annotated_genes_by_strand[chr][strand][elem])

#   for (sgid, name, gene_strand, gene_start, gene_end) in annotated_genes_by_strand[chr][strand]:
#        
#       if coordsOverlap(start, end, gene_start, gene_end):
#           if strand == ".":
#               gene_list.append(name)
#           elif gene_strand == strand:
#               gene_list.append(name)
#           elif strand is None:
#               gene_list.append(name)

#   if strand is None:
#       if len(gene_list) > 1:
#           return "None"
#       else:
#           return gene_list[0]
    gene_list = list(gene_set)

    return ",".join(gene_list)

def inferIE_junction(incl_junction, excl_junctions):
    """
    Returns the ie_junction coordinate.
    """
    chr, incl_start, incl_end = convertCoordStr(incl_junction)

    chr, excl_start, excl_end = convertCoordStr(getLongestJunction(excl_junctions.split(",")))

    if incl_start == excl_start:
        return formatCoordStr(chr, excl_end, excl_end + 1)
    
    # Else the ends are equal
    return formatCoordStr(chr, excl_start - 1, excl_start)

def inferInclusionJunctions(chr, excl_start, excl_end, incl_coord_list):
    coord_list = []

    for coord in incl_coord_list:
        chr, coord_start, coord_end = convertCoordStr(coord)

        coord_list.append((coord_start, coord_end))

    coord_list.sort()     

    prev_start = excl_start
    
    jcns = []
    for coord in coord_list:
        this_jcn = formatCoordStr(chr, prev_start, coord[0] - 1)
        jcns.append(this_jcn)
        prev_start = coord[1] + 1

    # Add last junction
    last_jcn = formatCoordStr(chr, prev_start, excl_end)

    jcns.append(last_jcn)

    return jcns

def isAnnotated(coord_str, annotated_coords):
    """ 
    coord_str - chr_start_end
    annotated_coords - {chr:set([(start, end, strand)])}
    """
    (chr, start, end) = convertCoordStr(coord_str)

    if chr in annotated_coords:
        if (start, end, "+") in annotated_coords[chr]:
            return True
        if (start, end, "-") in annotated_coords[chr]:
            return True

    return False

def isConnected(exon_list, chr, coord_start2end):
    """
    Checks if an intron connects any of the exons.
    """
    exon_list.sort()
    start_idx = 0
    end_idx = 1
    for i in range(len(exon_list)):
        for j in range(len(exon_list)):
            if i >= j:
                continue
      
            intron_start = exon_list[i][end_idx] + 1
            intron_end = exon_list[j][start_idx] - 1 
            # check for intron
            if intron_start in coord_start2end[chr]:
                if intron_end in coord_start2end[chr][intron_start]:
                    return True
            
    return False
            
def normalizeByLen(count, seq_len):
    if not DO_LEN_NORM:
        return count

    return int(round(count/(seq_len/DEF_EXON_LEN_NORM)))

def parse_all_as_event_regions(line_elem):

    returnList1 = line_elem.split(";")

    returnList2 = []
    for item1 in returnList1:
        if item1 == "":
            continue
        # Further separate by ","
        returnList2.extend(item1.split(","))

    return returnList2

def parse_qnamefile(region2qname_file_name):
    """
    Returns dictionary with information from the file.
    """    
    file = open(region2qname_file_name)
    
    region2qname2count = {}

    for line in file:
        line = formatLine(line)

        lineList = line.split("\t")

#        region_list = lineList[0].split()
        region_list = list(convertCoordStr(lineList[0]))
        region_list[0] = formatChr(region_list[0])

#        region = .join(region_list)
        region = "%s:%d-%d" % (region_list[0],
                               region_list[1],
                               region_list[2])

        if region in region2qname2count:
            # Debugging...this shouldn't happen
            ERROR_LOG.write("Region occurs multiple times in file: %s:%s\n" % (region2qname_file_name,
                                                                             region))
            continue

        qnameList = lineList[1].split(",")

        region2qname2count[region] = {}

        for qname in qnameList:
            if qname in region2qname2count[region]:
                region2qname2count[region][qname] += 1
            else:
                region2qname2count[region][qname] = 1

    file.close()

    return region2qname2count


def parseCoordCounts(mapped_file_name, norm):
    """
    Returns dictionary {coord: count}
    """
    coord2count = {}

    file = open(mapped_file_name)

    for line in file:
        line = formatLine(line)
        chr, start, end, count_str = line.split("\t")

        chr = formatChr(chr)

        k = formatCoordStr(chr, start, end)

        count = int(count_str)

        if norm:
            count = int(round(count/norm))

        coord2count[k] = count

    return coord2count

def parseIEJcnFiles(ie1_file, ie2_file):
    """ 
    Returns a dictionary {jcn_coord_str: {"left":[file1 count, file2 count],
                                          "right":[file1_count, file2_count]}
    """
    ie_jctn_dict = {}

    for line in ie1_file:
        line = formatLine(line)

        jcn_str, left_count_str, right_count_str = line.split("\t")

#        chr, start, end = jcn_str.split()
        (chr, start, end) = convertCoordStr(jcn_str)
        chr = formatChr(chr)
#        jcn_str = .join([chr, start,end])
        jcn_str = formatCoordStr(chr, start, end)
        
        left_count = int(left_count_str)
        right_count = int(right_count_str) 

        ie_jctn_dict[jcn_str] = {"left": [left_count, 0],
                                 "right":[right_count,0]}

    for line in ie2_file:
        line = formatLine(line)

        jcn_str, left_count_str, right_count_str = line.split("\t")

#        chr, start, end = jcn_str.split()
        chr, start, end = convertCoordStr(jcn_str)
        chr = formatChr(chr)
#        jcn_str = .join([chr,start,end])
        jcn_str = formatCoordStr(chr, start, end)

        left_count = int(left_count_str)
        right_count = int(right_count_str)

        if jcn_str in ie_jctn_dict:
            ie_jctn_dict[jcn_str]["left"][1] = left_count
            ie_jctn_dict[jcn_str]["right"][1] = right_count
        else:
            ie_jctn_dict[jcn_str] = {"left": [0, left_count],
                                     "right":[0, right_count]}

    return ie_jctn_dict

def parseJcns(jcn_file1, jcn_file2, genome_file, disambiguate_jcn_strand):
    jcn_count_dict = {}
    coord_start2end = {}
    coord_end2start = {}

    # {chr_start_end: strand}
    jcn2strand = {}

    # {chr: strand: coord search tree}
    jcn_search_tree = {}

    jcnFilesHaveStrandInfo = False

    for line in jcn_file1:
        line = formatLine(line)

        if line.startswith("track"):
            continue

        # Allowing junctions to == 0.  
        # Enforcing that counts have to be greater than 0
#        if int(count) <= 0:
#            continue

        # coord_str: chr_start_end
        coord_str, strand, count = translateInput(line)
    
        if coord_str is None:
            continue

        jcn2strand[coord_str] = strand

        # Initializing count dictionary
        jcn_count_dict[coord_str] = [int(count), 0]

        (chr, start, end) = convertCoordStr(coord_str)

        chr = formatChr(chr)

        # Setting up start2end
        if chr in coord_start2end:
            updateDictOfLists(coord_start2end[chr], int(start), int(end))
        else:
            coord_start2end[chr] = {int(start):[int(end)]}
        
        # Setting up end2start
        if chr in coord_end2start:
            updateDictOfLists(coord_end2start[chr], int(end), int(start))
        else:
            coord_end2start[chr] = {int(end):[int(start)]}

    # Now do file 2
    for line in jcn_file2:
        line = formatLine(line)

        if line.startswith("track"):
            continue

        # Allowing junction to be equal to 0
        # Enforcing that counts have to be greater than 0
#        if int(count) <= 0:
#            continue

        # coord_str: chr_start_end
        coord_str, strand, count = translateInput(line)

        if coord_str is None:
            continue

        try:
            jcn2strand[coord_str] = updateStrand(jcn2strand[coord_str], strand)
        except:
            jcn2strand[coord_str] = strand

        if coord_str in jcn_count_dict:
            jcn_count_dict[coord_str][1] = int(count)
        else:
            jcn_count_dict[coord_str] = [0, int(count)]
            
            (chr, start, end) = convertCoordStr(coord_str)

            chr = formatChr(chr)

            # Adding new coords. 
            if chr in coord_start2end:
                updateDictOfLists(coord_start2end[chr], int(start), int(end))
            else:
                coord_start2end[chr] = {int(start):[int(end)]}

            if chr in coord_end2start:
                updateDictOfLists(coord_end2start[chr], int(end), int(start))
            else:
                coord_end2start[chr] = {int(end):[int(start)]}

    # Populate jcn2strand dict for strands that are unknown
    if disambiguate_jcn_strand:
        getJcnStrandInfo(jcn2strand, genome_file)

    # Build jcn_search from this structure
    jcn_dict = {} # Prior to making search tree
    for jcn_str in jcn2strand:
        chr, start, end = convertCoordStr(jcn_str)
        
        if chr in jcn_dict:
            updateDictOfLists(jcn_dict[chr], jcn2strand[jcn_str], (start, end))
        else: 
            jcn_dict[chr] = {jcn2strand[jcn_str]: [(start,end)]}
    # Build search tree
    for chr in jcn_dict:
        jcn_search_tree[chr] = {}
        for strand in jcn_dict[chr]:
            jcn_search_tree[chr][strand] = getSearchTree(jcn_dict[chr][strand])

        # Add empty search trees for remaining strands
        for strand in ["+", "-", "."]:
            if strand not in jcn_search_tree[chr]:
                jcn_search_tree[chr][strand] = getSearchTree([])
    
    return (jcn_count_dict, coord_start2end, coord_end2start, jcn2strand, jcn_search_tree)

def parseReadAssocFile(paired_read_w_coord_file_name):
    """ 
    Returns a dictionary that maps the coordinate to a dictionary containing
    qnames to the number of occurences.
    {coord:{qname:count}}
    """
    file = open(paired_read_w_coord_file_name)

    # One read can have multiple counts to a region if they are paired end or
    # even strobe reads.
    coord2read2count = {}

    for line in file: 

        line = formatLine(line)
        lineList = line.split("\t")

#        coord_comp = lineList[-1].split()
        coord_comp = convertCoordStr(lineList[-1])

        coord_comp[0] = formatChr(coord_comp[0])

#        coord = .join(coord_comp)
        coord = formatCoordStr(coord_comp[0],
                               coord_comp[1],
                               coord_comp[2])

        qname = lineList[0]

        if coord in coord2read2count:
            if qname in coord2read2count[coord]:
                coord2read2count[coord][qname] += 1
            else: 
                coord2read2count[coord][qname] = 1
        else:
            coord2read2count[coord] = {qname:1}

    file.close()

    return coord2read2count

def printAlternativeDonorsAcceptors(db,
                                    annotated_genes,
                                    annotated_genes_by_strand,
                                    alt_first_exons,
                                    alt_last_exons,
                                    alt_first_exons_start2end,
                                    alt_first_exons_end2start,
                                    alt_last_exons_start2end,
                                    alt_last_exons_end2start,
                                    all_jcn_count_dict, 
                                    param_all_coord_start2end,
                                    param_all_coord_end2start, 
                                    all_jcn2strand,
                                    ir_count_dict,
                                    cassette_exons,
                                    annotated_introns,
                                    annotated_exons,
                                    annotated_exons_by_strand,
                                    all_confident_exons,
                                    all_confident_exons_start2end,
                                    all_confident_exons_end2start,
                                    printExonCoord,
                                    exon_coords,
                                    excl_jcns,
                                    donor_out, afe_out, jcn_only_donor_out,
                                    accept_out, ale_out, jcn_only_accept_out,
                                    all_event_info_out,
                                    norm1, norm2):

    (all_coord_start2end, 
     all_coord_end2start) = filterOutExclJcns(param_all_coord_start2end,
                                              param_all_coord_end2start,
                                              excl_jcns)

    filtered_cassette_exons_by_end = copyCassetteDict(cassette_exons)
    filtered_cassette_exons_by_start = copyCassetteDict(cassette_exons)
    removeOverlappingCassetteExons(filtered_cassette_exons_by_end, "end")
    removeOverlappingCassetteExons(filtered_cassette_exons_by_start, "start")

    filtered_cassette_exons_by_end_end2start = {}
    for chr in filtered_cassette_exons_by_end:
        for (start, end) in filtered_cassette_exons_by_end[chr]:
            updateDictExons(filtered_cassette_exons_by_end_end2start,
                            chr, end, start)
    filtered_cassette_exons_by_start_start2end = {}
    for chr in filtered_cassette_exons_by_start:
        for (start, end) in filtered_cassette_exons_by_start[chr]:
            updateDictExons(filtered_cassette_exons_by_start_start2end,
                            chr, start, end)
                

    for chr in all_coord_end2start:
        for end in all_coord_end2start[chr]:
            # Check if there are no alternative donors or acceptors
            if len(all_coord_end2start[chr][end]) == 1:
                continue

            # Make sure all strands agree
            strand = None
            onlyNovelIntrons = True
            for start in all_coord_end2start[chr][end]:
                this_jcn = formatCoordStr(chr, start, end)
                strand = updateStrand(strand, all_jcn2strand[this_jcn])
                if chr in annotated_introns:
                    if (start, end, "+") in annotated_introns[chr]:
                        onlyNovelIntrons = False
                    elif (start, end, "-") in annotated_introns[chr]:
                        onlyNovelIntrons = False
                    else: # A novel intron
                        continue

            if strand is "None":
                continue

            # Check if the upstream exon is a cassette exon
            if hasAdjExons(chr, filtered_cassette_exons_by_end_end2start, all_coord_end2start[chr][end],"P"):
                continue

            # Check if the event is simple, meaning it is not an alternative
            # first exon event
        
            # Becuase I made changes to AFE and ALE events where exons that
            # overlap are clustered into the same isoform, I had to add an
            # additional loop in the main part of this function
            ad_aa_afe_ale_events = [{"isSimple": True,
                                     "isAltFirstLast": False,
                                     "all_coord_end2start":{chr:
                                                            {end:
                                                             list(all_coord_end2start[chr][end])}}}]

            if strand == "+":
                if annotated_exons:
                    overlap_result = exonsOverlapping(chr, all_confident_exons,
                                                      all_confident_exons_start2end,
                                                      all_confident_exons_end2start,
                                            all_coord_end2start[chr][end],"P")
                    if overlap_result == False:
                        ad_aa_afe_ale_events[0]["isSimple"] = False
                        if hasAdjExons(chr, alt_first_exons_end2start, all_coord_end2start[chr][end],
                                       "P"):
                            event_jcns = []
                            for this_s in all_coord_end2start[chr][end]:
                                event_jcns.append(formatCoordStr(chr, this_s,
                                                                end))

                            # ad_aa_afe_ale_events list will be updated
                            find_AFE_ALE_clusters(ad_aa_afe_ale_events,
                                                  "P",
                                                  all_confident_exons,
                                                  all_confident_exons_start2end,
                                                  all_confident_exons_end2start,
                                                  all_jcn_count_dict,
                                                  chr, 
                                                  event_jcns)
                    elif overlap_result is None:
                        ad_aa_afe_ale_events[0]["isSimple"] = False
        
            else:
                if annotated_exons:
                    overlap_result = exonsOverlapping(chr, all_confident_exons,
                                                      all_confident_exons_start2end,
                                                      all_confident_exons_end2start,
                                                      all_coord_end2start[chr][end],
                                                      "P")
                    if overlap_result == False:
                        ad_aa_afe_ale_events[0]["isSimple"] = False
                        if hasAdjExons(chr, alt_last_exons_end2start, all_coord_end2start[chr][end],
                                       "P"):
                            event_jcns = []
                            for this_s in all_coord_end2start[chr][end]:
                                event_jcns.append(formatCoordStr(chr, this_s,
                                                                end))

                            # ad_aa_afe_ale_events list will be updated
                            find_AFE_ALE_clusters(ad_aa_afe_ale_events,
                                                  "P",
                                                  all_confident_exons,
                                                  all_confident_exons_start2end,
                                                  all_confident_exons_end2start,
                                                  all_jcn_count_dict,
                                                  chr, 
                                                  event_jcns)

                    elif overlap_result is None:
                        ad_aa_afe_ale_events[0]["isSimple"] = False

            # BEGIN OF MAIN LOOP
            for event_dict in ad_aa_afe_ale_events:
                isSimple = event_dict["isSimple"]
                isAltFirstLast = event_dict["isAltFirstLast"]
            
                this_all_coord_end2start = event_dict["all_coord_end2start"]
            
                # Check if any introns are novel and which position is
                # longest
                isNovel = False
                longest_start = INFINITY
                for start in this_all_coord_end2start[chr][end]:
                    if start < longest_start:
                        longest_start = start
                   
                    this_distal_jcn = formatCoordStr(chr,
                                                  start, end) 
                    # Check if it is a novel intron
                    if chr not in annotated_introns:
                        isNovel = True
                    else:
                        if isAltFirstLast:
                            for jcn_str in event_dict["jcn2jcn_str"][this_distal_jcn].split(","):
                                c, j_start, j_end = convertCoordStr(jcn_str)
                                if (j_start, j_end, strand) not in annotated_introns[chr]:
                                    isNovel = True
                        else:
                            if (start, end, strand) not in annotated_introns[chr]:
                                isNovel = True


                # Make parallel count dict {start: (count_raw, count_lenNorm)}
                par_jcn_count_dict = {}
                for start in this_all_coord_end2start[chr][end]:
                    jcn_coord_str = formatCoordStr(chr, start, end)
                    if isAltFirstLast:
                        par_jcn_count_dict[start] = event_dict["jcn_cluster_sum"][jcn_coord_str]
                    else:
                        par_jcn_count_dict[start] = all_jcn_count_dict[jcn_coord_str]

                # Get counts and output
                for start in this_all_coord_end2start[chr][end]:
                    # If only two isoforms, then only calculation not performed for
                    # longest intron
                    if len(this_all_coord_end2start[chr][end]) == 2:
                        if start == longest_start:
                            continue
                
                    this_distal_jcn = formatCoordStr(chr, start, end)

                    exclusion_raw = 0
                    inclusion_raw = 0
#                   exclusion_lenNorm = 0
#                   inclusion_lenNorm = 0

                    inclusion_jcn_raw = 0
#                    inclusion_jcn_lenNorm = 0

                    exclusion_str_list = []
                    exclusion_distal_list = []
                    exclusion_jcn_raw_list = []
#                    exclusion_jcn_lenNorm_list = []

                    # Add novel junctions to Alt First and last exons
                    if isAltFirstLast:
                        if event_dict["novel_jcns"]:
                            exclusion_str_list.append(event_dict["novel_jcns"])
                            
                            exclusion_distal_list.append(event_dict["novel_jcns"].split(",")[0])

                            novel_sum_raw = event_dict["novel_jcn_sum_samp2"]

                            if norm2:
                                novel_sum_raw = int(round(novel_sum_raw/norm2))

                            exclusion_jcn_raw_list.append(novel_sum_raw)

                    exclusion_exon_raw = 0
                    inclusion_exon_raw = 0
#                   exclusion_exon_cts2 = 0
#                   inclusion_exon_cts2 = 0
               
                    ie_jcn = None 
                    const_regions = []
                    const_str = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, end + 1, "N")
                    if const_str:
                        const_regions.append(const_str)
                        if printExonCoord:
                            exon_coords.add(convertCoordStr(const_str))
                   
                    for par_start in par_jcn_count_dict:
                        if par_start == start:
    #                       # No paired end counting yet
#                            this_jcn_raw = par_jcn_count_dict[par_start][0]
                            this_jcn_raw = par_jcn_count_dict[par_start][1]

                            if norm2:
                                this_jcn_raw = int(round(this_jcn_raw/norm2))

                            # Not normalizing by length here because the
                            # junction counts need to be maintained for
                            # downstream proportional analysis

                            inclusion_raw += this_jcn_raw
#                            inclusion_cts2 += this_jcn_lenNorm
            
                            inclusion_jcn_raw = this_jcn_raw
#                            inclusion_jcn_ct2 = this_jcn_lenNorm

                        else:
                            excl_intron = formatCoordStr(chr,par_start,end)
                            if isAltFirstLast:
                                exclusion_str_list.append(event_dict["jcn2jcn_str"][excl_intron])
                            else:
                                exclusion_str_list.append(excl_intron)

                            exclusion_distal_list.append(excl_intron)

#                            this_excl_ct1 = par_jcn_count_dict[par_start][0]
                            this_excl_raw = par_jcn_count_dict[par_start][1]

                            if norm2:
#                                this_excl_ct1 = int(round(this_excl_ct1/norm1))
                                this_excl_raw = int(round(this_excl_raw/norm2))

                            # Not length normalizing here to maintain junction
                            # proportions for later analysis
                            
                            exclusion_raw += this_excl_raw
#                            exclusion_cts2 += this_excl_ct2

                            exclusion_jcn_raw_list.append(this_excl_raw)
#                            exclusion_jcn_cts2_list.append(this_excl_ct2)


                    n_or_k = "K"
                    if isNovel:
                        n_or_k = "N"
                    

                    # ordered_pos is an ordered list.  The proportions variables are in the
                    # same order as the ordered list.
                    (ordered_pos, 
                     proportions1, 
                     proportions2,
                     not_used1,
                     not_used2) = getSSOrderAndProportions("alt_start",
                                                            start, end,
                                                            exclusion_distal_list,
                                                            exclusion_jcn_raw_list,
                                                            inclusion_raw,
                                                            exclusion_jcn_raw_list,
                                                            inclusion_raw)

                    ie_jcns = []
                    ie_jcns_raw = []
#                    ie_jcn_cts2 = []
                    # Add any intron retention counts associated with the
                    # intron/exon boundaries of the isoforms.
                    # If there are paired-end samples, counting does not occur
                    # until later in updateCounts2AltDonorAccept
                    if isSimple:
                        # Get all IE junction counts that correspond to all
                        # positions except the last
                        for i in range(len(ordered_pos)-1):
                            intron = formatCoordStr(chr, ordered_pos[i], end)
               
                            ie_jcn = formatCoordStr(chr, ordered_pos[i] - 1,
                                                   ordered_pos[i]) 
                            ie_jcns.append(ie_jcn)
                          
                            ie_jcn_raw = 0
#                            ie_jcn_ct2 = 0

                            if ir_count_dict:
                                if intron in ir_count_dict:

                                    ie_jcn_raw = ir_count_dict[intron]["left"][1]

                                    if norm2:
                                        ie_jcn_raw = int(round(ie_ct_raw/norm2))

                                    # length normalizing at later step

                                    ie_jcns_raw.append(ie_jcn_raw)
#                                    ie_jcn_cts2.append(left_ct_lenNorm)
                                else:
                                    ie_jcns_raw.append(0)

                    # Add Exclusion or Inclusion Annotation
                    # Checks percent exclusion from both files.
#                   e_or_i = checkExclusionInclusion_AA_AD_AFE_ALE("alt_start",
#                                                                  start, end,
#                                                                  ordered_pos,
#                                                                  proportions1,
#                                                                  proportions2)
                    gene_name = inferGeneName(annotated_genes_by_strand, chr, strand, 
                                              exclusion_str_list + [formatCoordStr(chr,start,end)])

                    out_str = "%s\t%s\t%s\t%s\t%s" % (n_or_k,
                                                  "?",
                                                  gene_name,
                                                  chr,
                                                  strand)

                    if isAltFirstLast:
                        out_str += "\t%s" % event_dict["jcn2jcn_str"][this_distal_jcn]
                    else:
                        out_str += "\t%d\t%d" % (start, end)
                
                    out_str += "\t%s\t%s\t%d\t%s\t%d" % (";".join(exclusion_str_list),
                                                         ";".join(map(repr,exclusion_jcn_raw_list)),
                                                         inclusion_raw,
                                                         ";".join(map(repr,exclusion_jcn_raw_list)),
                                                         inclusion_raw)
                    inclusion_region = None
                    excl_regions = []

                    if isSimple: 
                        inclusion_regions = breakInclusionRegion("alt_start",
                                                                 chr,
                                                                 ordered_pos)

                        inclusion_region = ";".join(inclusion_regions)

                        if printExonCoord:
                            out_str += "\t%s" % inclusion_region
                            for incl_reg in inclusion_regions:
                                exon_coords.add(convertCoordStr(incl_reg))


                        # Add an additional constitutive region
                        const_str = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, longest_start - 1, "P")
                        if const_str:
                            const_regions.append(const_str)
                            if printExonCoord:
                                exon_coords.add(convertCoordStr(const_str))
                    elif isAltFirstLast:
                        # INCLUSION REGION
                        exon_coord_str = event_dict["jcn2exon_str"][this_distal_jcn]

                        inclusion_region = exon_coord_str

                        if exon_coord_str and exon_coord_str != "None":
                            if printExonCoord:
                                out_str += "\t%s" % exon_coord_str
                                exon_coords.add(convertCoordStr(exon_coord_str))
                        else:
                            if printExonCoord:
                                out_str += "\tNone"

                        # EXCLUSION REGION
                        for excl_start in par_jcn_count_dict:
                            if excl_start == start:
                                continue

                            this_excl_jcn = formatCoordStr(chr, excl_start, end)
                            exon_coord_str = event_dict["jcn2exon_str"][this_excl_jcn]

                            if exon_coord_str and exon_coord_str != "None":
                                excl_regions.append(exon_coord_str)
                                if printExonCoord:
                                    exon_coords.add(convertCoordStr(exon_coord_str))

                        if printExonCoord:   
                            out_str += "\t" 
                            if excl_regions == []:
                                out_str += "None"
                            else:
                                out_str += ",".join(excl_regions)                                

                    else: # Is nonSimple and not an AFE/ALE.  Only use junction
                          # counts
                        if printExonCoord:
                            out_str += "\tNone"
                        
                    out_str += "\n"
                    type = ""
                    if strand == "+":
                        if isSimple:
                            donor_out.write(out_str)
                            jcn_only_donor_out.write(out_str)
                            type = "alternative_donor"
                        elif isAltFirstLast:
                            afe_out.write(out_str)
                            type = "alternative_first_exon"
                        else:
                            jcn_only_donor_out.write(out_str)
                            type = "jcn_only_AD"
                    else:
                        if isSimple:
                            accept_out.write(out_str)
                            jcn_only_accept_out.write(out_str)
                            type = "alternative_acceptor"
                        elif isAltFirstLast:
                            ale_out.write(out_str)
                            type = "alternative_last_exon"
                        else:
                            jcn_only_accept_out.write(out_str)
                            type = "jcn_only_AA"
        

                    if not inclusion_region:
                        inclusion_region = ""
                    if not ie_jcn:
                        ie_jcn = ""
                    
                    if isAltFirstLast:
                        inclusion_str = event_dict["jcn2jcn_str"][this_distal_jcn]
                    else:
                        inclusion_str = this_distal_jcn

                    const_str = ";".join(const_regions)

                    ie_jcns_raw_sum = 0
                    if len(ie_jcns_raw) >= 1:
                        ie_jcns_raw_sum = sum(ie_jcns_raw) 
    
                    out_str= getAllEventStr(n_or_k, type, gene_name, chr, strand, 
                                             ";".join(exclusion_str_list),
                                             inclusion_str,
                                             ";".join(excl_regions), inclusion_region,
                                             ";".join(ie_jcns),
                                             const_str,
                                             ";".join(map(repr,exclusion_jcn_raw_list)),
                                             repr(inclusion_jcn_raw),
                                             repr(exclusion_raw),
                                             repr(inclusion_jcn_raw),
                                             "","",
                                             None, None,
                                             ";".join(map(repr, ie_jcns_raw)),
                                             ie_jcns_raw_sum,
                                             "",
                                             None)
                    all_event_info_out.write(out_str + "\n")

                    # If the event is an alternative donor or acceptor type, then
                    # it also needs to be printed again as a jcn_only_* event
                    if type == "alternative_donor":
                        all_event_info_out.write(out_str.replace("alternative_donor",
                                                                 "jcn_only_AD") 
                                                 + "\n")
                    elif type == "alternative_acceptor":
                        all_event_info_out.write(out_str.replace("alternative_acceptor",
                                                                 "jcn_only_AA")
                                                 + "\n")
                                           
    # Now look at the other dictionary            
    for chr in all_coord_start2end:
        for start in all_coord_start2end[chr]:
            # Check if there are no alternative donors or acceptors
            if len(all_coord_start2end[chr][start]) == 1:
                continue

            # Make sure all strands agree
            strand = None
            onlyNovelIntrons = True
            for end in all_coord_start2end[chr][start]:
                this_jcn = formatCoordStr(chr, start, end)
                strand = updateStrand(strand, all_jcn2strand[this_jcn])
                if chr in annotated_introns:
                    if (start, end, "+") in annotated_introns[chr]:
                        onlyNovelIntrons = False
                    elif (start, end, "-") in annotated_introns[chr]:
                        onlyNovelIntrons = False
                    else: # A novel intron
                        continue

            if strand is "None":
                continue

            # Check if the downstream exon is a cassette exon
            if hasAdjExons(chr, filtered_cassette_exons_by_start_start2end,
                              all_coord_start2end[chr][start],"N"):
                continue

            ad_aa_afe_ale_events = [{"isSimple": True,
                                     "isAltFirstLast": False,
                                     "all_coord_start2end":{chr:
                                                            {start:
                                                             list(all_coord_start2end[chr][start])}}}]

            # Check if the event is simple, meaning it is not an alternative
            # first or last exon event
            if strand == "-":
                if annotated_exons:
                    overlap_result = exonsOverlapping(chr, all_confident_exons,
                                                      all_confident_exons_start2end,
                                                      all_confident_exons_end2start,
                                                      all_coord_start2end[chr][start],
                                                      "N")
                    if overlap_result == False:
                        ad_aa_afe_ale_events[0]["isSimple"] = False
                        if hasAdjExons(chr, alt_first_exons_start2end,
                                       all_coord_start2end[chr][start], "N"):
                            event_jcns = []
                            for this_e in all_coord_start2end[chr][start]:
                                event_jcns.append(formatCoordStr(chr, start,
                                                                this_e))
        
                            # ad_aa_afe_ale_events list will be updated
                            find_AFE_ALE_clusters(ad_aa_afe_ale_events,
                                                  "N",
                                                  all_confident_exons,
                                                  all_confident_exons_start2end,
                                                  all_confident_exons_end2start,
                                                  all_jcn_count_dict,
                                                  chr,
                                                  event_jcns)
                    elif overlap_result is None:
                        ad_aa_afe_ale_events[0]["isSimple"] = False

            else: # strand is +
                if annotated_exons:
                    overlap_result = exonsOverlapping(chr, all_confident_exons,
                                                      all_confident_exons_start2end,
                                                      all_confident_exons_end2start,
                                                      all_coord_start2end[chr][start],
                                                      "N")
                    if overlap_result == False:
                        ad_aa_afe_ale_events[0]["isSimple"] = False
                        if hasAdjExons(chr, alt_last_exons_start2end,
                                       all_coord_start2end[chr][start],
                                       "N"):
                            event_jcns = []
                            for this_e in all_coord_start2end[chr][start]:
                                event_jcns.append(formatCoordStr(chr, start,
                                                                this_e))

                            # ad_aa_afe_ale_events list will be updated
                            find_AFE_ALE_clusters(ad_aa_afe_ale_events,
                                                  "N",
                                                  all_confident_exons,
                                                  all_confident_exons_start2end,
                                                  all_confident_exons_end2start,
                                                  all_jcn_count_dict,
                                                  chr, 
                                                  event_jcns)
                                
                    elif overlap_result is None:
                        ad_aa_afe_ale_events[0]["isSimple"] = False

            # BEGIN OF MAIN LOOP
            for event_dict in ad_aa_afe_ale_events:    
                isSimple = event_dict["isSimple"]
                isAltFirstLast = event_dict["isAltFirstLast"]

                this_all_coord_start2end = event_dict["all_coord_start2end"]        

                # Check if any introns are novel and which position is
                # longest...thus always treated as the "exclusion" event
                isNovel = False
                longest_end = -1
                for end in this_all_coord_start2end[chr][start]:
                    if end > longest_end:
                        longest_end = end

                    this_distal_jcn = formatCoordStr(chr,start, end)
                    
                    # Check if it is a novel intron
                    if chr not in annotated_introns:
                        isNovel = True
                    else:
                        if isAltFirstLast:
                            for jcn_str in event_dict["jcn2jcn_str"][this_distal_jcn].split(","):
                                c, j_start, j_end = convertCoordStr(jcn_str)
                                if (j_start, j_end, strand) not in annotated_introns[chr]:
                                    isNovel = True
                        else:
                            if (start, end, strand) not in annotated_introns[chr]:
                                isNovel = True

                # Make parallel count dict {end: (count_raw, count_lenNorm)}
                par_jcn_count_dict = {}
                for end in this_all_coord_start2end[chr][start]:
                    jcn_coord_str = formatCoordStr(chr, start, end)
                    if isAltFirstLast:
                        par_jcn_count_dict[end] = event_dict["jcn_cluster_sum"][jcn_coord_str]
                    else:
                        par_jcn_count_dict[end] = all_jcn_count_dict[jcn_coord_str]

                # Get counts and output
                for end in this_all_coord_start2end[chr][start]:
                    # If there are only two isoforms, then calculation not
                    # performed for longest intron
                    if len(this_all_coord_start2end[chr][start]) == 2:
                        if end == longest_end:
                            continue

                    this_distal_jcn = formatCoordStr(chr, start, end)

                    exclusion_raw = 0
                    inclusion_raw = 0

                    inclusion_jcn_raw = 0

                    exclusion_str_list = []
                    exclusion_distal_list = []
                    exclusion_jcn_raw_list = []

                    # Add novel junctions to alt first and last exons
                    if isAltFirstLast:
                        if event_dict["novel_jcns"]:
                            exclusion_str_list.append(event_dict["novel_jcns"])

                            exclusion_distal_list.append(event_dict["novel_jcns"].split(",")[0])
        
                            novel_sum_raw = event_dict["novel_jcn_sum_samp2"]

                            if norm2:
                                novel_sum_raw = int(round(novel_sum_raw/norm2))

                            exclusion_jcn_raw_list.append(novel_sum_raw)

                    exclusion_exon_raw = 0
                    inclusion_exon_raw = 0
       
                    ie_jcn = None 
                    const_regions = []
                    const_str = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, start - 1, "P")
                    if const_str:
                        const_regions.append(const_str)
                        if printExonCoord:
                            exon_coords.add(convertCoordStr(const_str))

                    for par_end in par_jcn_count_dict:
                        if par_end == end:
                            # No paired end counting yet
                            this_incl_raw = par_jcn_count_dict[par_end][1]

                            if norm2:
                                this_incl_raw = int(round(this_incl_raw/norm2))

                            # Not length normalizing to maintain proportions
                            # for later calculation
        
                            inclusion_raw += this_incl_raw

                            inclusion_jcn_raw = this_incl_raw

                        else:
                            excl_intron = formatCoordStr(chr, start, par_end)
    
                            if isAltFirstLast:
                                exclusion_str_list.append(event_dict["jcn2jcn_str"][excl_intron])
                            else:
                                exclusion_str_list.append(excl_intron)
        
                            exclusion_distal_list.append(excl_intron)

                            this_raw = par_jcn_count_dict[par_end][1]

                            if norm2:
                                this_raw = int(round(this_raw/norm2))

                            # Not length normalizing here to maintain junction
                            # proportions for later step

                            exclusion_raw += this_raw
                            
                            exclusion_jcn_raw_list.append(this_raw)

                    n_or_k = "K"
                    if isNovel:
                        n_or_k = "N"

                    # ordered_pos is an ordered list.  The proportions variables are in the
                    # same order as the ordered list.
                    (ordered_pos, 
                     proportions1, 
                     proportions2,
                     not_used1,
                     not_used2) = getSSOrderAndProportions("alt_end",
                                                           start, end,
                                                           exclusion_distal_list,
                                                           exclusion_jcn_raw_list,
                                                           inclusion_raw,
                                                           exclusion_jcn_raw_list,
                                                           inclusion_raw)
                    
                    # Add any intron retention counts associated with this
                    # intron/exon boundaries of the isoforms
                    # IE counts are added later if the reads are paired-end
                    ie_jcns = []
                    ie_jcns_raw = []
                    if isSimple:
                        # Get all IE junctions that correspond to all positions
                        # except the first
                        for i in range(1,len(ordered_pos)):
                            intron = formatCoordStr(chr, start, ordered_pos[i])

                            ie_jcn = formatCoordStr(chr, ordered_pos[i],
                                                   ordered_pos[i] + 1)
                            ie_jcns.append(ie_jcn)

                            ie_jcn_raw = 0

                            if ir_count_dict:
                                if intron in ir_count_dict:
                                    ie_jcn_raw = ir_count_dict[intron]["right"][1]

                                    if norm2:
                                        ie_jcn_raw = int(round(ie_jcn_raw/norm2))

                                    # not length normalizing here, will perform
                                    # normalization later
                                    ie_jcns_raw.append(ie_jcn_raw)
                                else:
                                    ie_jcns_raw.append(0)

                    # Add Exclusion or Inclusion Annotation
#                   e_or_i = checkExclusionInclusion_AA_AD_AFE_ALE("alt_end",
#                                                                  start, end,
#                                                                  ordered_pos,
#                                                                  proportions1,
#                                                                  proportions2)
                    gene_name = inferGeneName(annotated_genes_by_strand, chr, strand,
                                              exclusion_str_list + [inclusion_str]) 


                    out_str = "%s\t%s\t%s\t%s\t%s" % (n_or_k,
                                                      "?",
                                                      gene_name,
                                                      chr,
                                                      strand)
                                                     
                    if isAltFirstLast:
                        out_str += "\t%s" % event_dict["jcn2jcn_str"][this_distal_jcn]                       
                    else:
                        out_str += "\t%d\t%d" % (start, end)

                    out_str += "\t%s\t%s\t%d\t%s\t%d" % (";".join(exclusion_str_list),
                                                         ";".join(map(repr,exclusion_jcn_raw_list)),
                                                         inclusion_raw,
                                                         ";".join(map(repr,exclusion_jcn_raw_list)),
                                                         inclusion_raw)

                    inclusion_region = None
                    excl_regions = []

                    if isSimple:
                        inclusion_regions = breakInclusionRegion("alt_end",
                                                                 chr,
                                                                 ordered_pos)

                        inclusion_region = ";".join(inclusion_regions)

                        if printExonCoord:
                            out_str += "\t%s" % inclusion_region
               
                            for incl_reg in inclusion_regions: 
                                exon_coords.add(convertCoordStr(incl_reg))

                        # Add an additional constituve region
                        const_str = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, longest_end + 1, "N")
                        if const_str:
                            const_regions.append(const_str)
                            if printExonCoord:
                                exon_coords.add(convertCoordStr(const_str))

                    elif isAltFirstLast:
                        exon_coord_str = event_dict["jcn2exon_str"][this_distal_jcn] 

                        inclusion_region = exon_coord_str

                        if exon_coord_str and exon_coord_str != "None":
                            if printExonCoord:
                                out_str += "\t%s" % exon_coord_str
                                exon_coords.add(convertCoordStr(exon_coord_str))
                        else:
                            if printExonCoord:
                                out_str += "\tNone"

                        # EXCLUSION REGION
                        for excl_end in par_jcn_count_dict:
                            if excl_end == end:
                                continue

                            this_excl_jcn = formatCoordStr(chr, start, excl_end)
                            exon_coord_str = event_dict["jcn2exon_str"][this_excl_jcn]

                            if exon_coord_str and exon_coord_str != "None":
                                excl_regions.append(exon_coord_str)
                                if printExonCoord:
                                    exon_coords.add(convertCoordStr(exon_coord_str))

                        if printExonCoord:
                            out_str += "\t"
                            if excl_regions == []:
                                out_str += "None"
                            else:
                                out_str += ",".join(excl_regions)

                    else: # is a jcn_only event
                        if printExonCoord:
                            out_str += "\tNone"

                    out_str += "\n"

                    type = ""
                    if strand == "-":
                        if isSimple:
                            donor_out.write(out_str)
                            jcn_only_donor_out.write(out_str)
                            type = "alternative_donor"
                        elif isAltFirstLast:
                            afe_out.write(out_str)
                            type = "alternative_first_exon"
                        else:
                            jcn_only_donor_out.write(out_str)
                            type = "jcn_only_AD"
                    else: # strand == +
                        if isSimple:
                            accept_out.write(out_str)
                            jcn_only_accept_out.write(out_str)
                            type = "alternative_acceptor"
                        elif isAltFirstLast:
                            ale_out.write(out_str)
                            type = "alternative_last_exon"
                        else:
                            jcn_only_accept_out.write(out_str)
                            type = "jcn_only_AA"


                    if not inclusion_region:
                        inclusion_region = ""
                    if not ie_jcn:
                        ie_jcn = ""

                    if isAltFirstLast:
                        inclusion_str = event_dict["jcn2jcn_str"][this_distal_jcn]
                    else:
                        inclusion_str = this_distal_jcn

                    const_str = ";".join(const_regions)

                    ie_jcns_raw_sum = 0
                    if len(ie_jcns_raw) >= 1:
                        ie_jcns_raw_sum = sum(ie_jcns_raw) 

                    out_str = getAllEventStr(n_or_k, type, gene_name, chr, strand, 
                                             ";".join(exclusion_str_list),
                                             inclusion_str, 
                                             ";".join(excl_regions), inclusion_region,
                                             ";".join(ie_jcns),
                                             const_str,
                                             ";".join(map(repr,exclusion_jcn_raw_list)),
                                             repr(inclusion_jcn_raw),
                                             repr(exclusion_raw),
                                             repr(inclusion_jcn_raw),
                                             "","",
                                             None,None,
                                             ";".join(map(repr, ie_jcns_raw)),
                                             ie_jcns_raw_sum, 
                                             "",
                                             None)
                    all_event_info_out.write(out_str + "\n")

                    if type == "alternative_donor":
                        all_event_info_out.write(out_str.replace("alternative_donor",
                                                                 "jcn_only_AD") 
                                                 + "\n")
                    elif type == "alternative_acceptor":
                        all_event_info_out.write(out_str.replace("alternative_acceptor",
                                                                 "jcn_only_AA")
                                                 + "\n")

def printAlternativePolyA(db, txt_db, 
                          annotated_genes,
                          full_exon_count_dict, 
                          full_multi_exon_count_dict, 
                          printExonCoords,
                          exon_coords,
                          read_length,
                          polya_out,
                          norm1, norm2):
    """
    Uses the alternative splicing database to identify the alternative polyA
    sites and then uses the annotations to count the exonic reads.
    """
    select = """SELECT DISTINCT splicing_graph2chr.splicing_graph_id,
                                splicing_graph2chr.chr,
                                splicing_graph2chr.strand,
                                exon.start,
                                exon.end
                FROM exon, node, splicing_graph2chr,
                     alternative_transcription_termination_site
                WHERE alternative_transcription_termination_site.node_id=node.node_id
                      AND
                      node.splicing_graph_id=splicing_graph2chr.splicing_graph_id
                      AND
                      node.node_id=exon.node_id"""

    records = db.getDBRecords_Dict(select, txt_db)

    # dictionary of putative alternative polyAs
    # { (chr, sgid, <anchor (start, or end depending on strand>): [exon_str,]
    alt_polyA_dict = {} 
    for row in records:
        sgid = int(row["splicing_graph_id"])
        chr = "chr" + row["chr"]
        start = int(row["start"])
        end = int(row["end"])
        strand = row["strand"]

        # Skip alternative last exons which involve alternative splicing
#        if chr in alt_last_exons:
#            if (start, end) in alt_last_exons[chr]:
#                continue

        if strand == "+":
            updateDictOfLists(alt_polyA_dict, (chr, sgid, start), 
                              formatCoordStr(chr, start, end))
        else:
            updateDictOfLists(alt_polyA_dict, (chr, sgid, end), 
                              formatCoordStr(chr, start, end))

    # Now find polyA events and print
    for (chr,sgid,anchor) in alt_polyA_dict:
        if len(alt_polyA_dict[(chr,sgid,anchor)]) < 2:
            continue

        # Prints difference between alt polyA exons
        diff, shortest_exon_str = getExonDistanceDifference(alt_polyA_dict[(chr,
                                                                            sgid,
                                                                            anchor)])
        # Only output alternative polyA events that have the potential to
        # notice a difference.
        if diff <= read_length:
            continue

        if hasOtherOverlappingGene(annotated_genes, alt_polyA_dict, chr, sgid, anchor):
            continue

        # Print out for each event
        for i in range(len(alt_polyA_dict[(chr,sgid,anchor)])): 
            # Only want to print out for inclusion events
            if alt_polyA_dict[(chr,sgid,anchor)][i] == shortest_exon_str:
                continue

            incl_file1_count = 0
            incl_file2_count = 0
            excl_file1_count = 0
            excl_file2_count = 0

            exon_str = alt_polyA_dict[(chr,sgid,anchor)][i]
           
            if full_exon_count_dict: # check if the dictionary was given
                # Ony read unique the exon are in the full_exon_count dict
                if exon_str in full_exon_count_dict: 
                    incl_file1_count = full_exon_count_dict[exon_str][0]
                    incl_file2_count = full_exon_count_dict[exon_str][1]
#                    excl_file1_count = full_exon_count_dict[exon_str][0]
#                    excl_file2_count = full_exon_count_dict[exon_str][1]
                # Reads that hit the exon, but not uniquely will be in the
                # multi_exon_count_dict.  These counts will be added the
                # exclusion file only
                if exon_str in full_multi_exon_count_dict:
                    excl_file1_count += full_multi_exon_count_dict[exon_str][0]
                    excl_file2_count += full_multi_exon_count_dict[exon_str][1]

        
            # Add Exclusion or Inclusion Annotation
            # Checks percent exclusion from both files.
            e_or_i = checkExclusionInclusion(excl_file1_count,
                                            incl_file1_count,
                                            excl_file2_count,
                                            incl_file2_count) 

           
#            exon_chr, exon_start, exon_end = convertCoordStr(exon_str.split(
            exon_chr, exon_start, exon_end = convertCoordStr(exon_str)
            gene_name = inferGeneName(annotated_genes_by_strand, chr, exon_start,
                        exon_end, ".")


            if norm1:
                excl_file1_count = int(round(excl_file1_count/norm1))
                incl_file1_count = int(round(incl_file1_count/norm1))
            
                excl_file2_count = int(round(excl_file2_count/norm2))
                incl_file2_count = int(round(incl_file2_count/norm2))
    

            # print
            out_str = "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d" % (e_or_i,
                                                            gene_name,
                                                            chr,
                                                            exon_str,
                                                            ",".join(alt_polyA_dict[(chr,sgid,anchor)]),
                                                            diff,
                                                            excl_file1_count,
                                                            incl_file1_count,
                                                            excl_file2_count,
                                                            incl_file2_count)                                                
            if printExonCoords:
                (inclusion_start, 
                 inclusion_end) = getInclusionPortion(exon_str,
                                                      alt_polyA_dict[(chr,sgid,anchor)])

                inclusion_str = formatCoordStr(chr, inclusion_start,
                                              inclusion_end)
                out_str += "\t%s\t" % inclusion_str
                exon_coords.add((chr, inclusion_start, inclusion_end))

                if inclusion_start > anchor:
                    exclusion_str = formatCoordStr(chr, 
                                                  anchor,
                                                  inclusion_start - 1)
                    exon_coords.add((chr, anchor, inclusion_start - 1))
                else:
                    exclusion_str = formatCoordStr(chr,
                                                  inclusion_end + 1,
                                                  anchor)
                    exon_coords.add((chr, inclusion_end + 1, anchor))
                out_str += exclusion_str
               
            polya_out.write(out_str + "\n")
            
def printCassetteExons(db,
                       alt_first_exons,
                       alt_last_exons,
                       annotated_genes,
                       annotated_genes_by_strand,
                       all_jcn_count_dict, 
                       all_coord_start2end, 
                       all_coord_end2start,
                       all_jcn2strand,
                       full_exon_count_dict,
                       full_multi_exon_count_dict,
                       start_multi_exon_count_dict,
                       end_multi_exon_count_dict,
                       annotated_exons,
                       annotated_exons_by_strand,
                       annotated_introns,
                       printExonCoord,
                       exon_coords,
                       excl_jcns,
                       cassette_out,
                       all_event_info_out,
                       norm1, norm2,
                       jcn_seq_len):
    """
    Checks for cassette exons by looking at exclusion, inclusion
    pairings.  Will note of the cassette exon is already annotated.
    """

    # {chr_exonstart_exonend:{"excl_set":set([chr_start_end]),
    #                          "left_set":set([chr_start_end]),
    #                          "right_set":set([chr_start_end])}
    cassette_exon_dict = {}
    
    return_cassette_exons = {}

    for chr in all_coord_start2end:
        for start in all_coord_start2end[chr]:
            if len(all_coord_start2end[chr][start]) > 1:

                all_coord_start2end[chr][start].sort()

                possible_excl = list(all_coord_start2end[chr][start][1:])

                strand = None

                for end in possible_excl:
                    if len(all_coord_end2start[chr][end]) > 1:
                        possible_left_ends = getPossibleLeftEnds(all_coord_start2end[chr][start],
                                                                 end)
                        for greater_start in all_coord_end2start[chr][end]:
                            for left_end in possible_left_ends:
                                if greater_start > left_end:
                                    excl_str = formatCoordStr(chr, start, end)
                                    strand = updateStrand(strand,
                                                          all_jcn2strand[excl_str])
                                    left_str = formatCoordStr(chr, start, left_end)
                                    strand = updateStrand(strand,
                                                          all_jcn2strand[left_str])
                                    right_str = formatCoordStr(chr,
                                                               greater_start,
                                                               end)
                                    strand = updateStrand(strand,
                                                          all_jcn2strand[right_str])

                                    exon_start = int(left_end) + 1
                                    exon_end = int(greater_start) - 1
                                    exon_coord = (chr, 
                                                  exon_start,
                                                  exon_end)

                                    # Enforce maximum length
                                    if exon_end - exon_start + 1 > MAX_EXON_LEN:
                                        continue

                                    # Make sure the exon is annotated
                                    if annotated_exons:
                                        exonFound = False
                                        if chr in annotated_exons:
                                            if (exon_start, exon_end, "+") in annotated_exons[chr]:
                                                exonFound = True
                                            if (exon_start, exon_end, "-") in annotated_exons[chr]:
                                                exonFound = True

                                        if not exonFound:
                                            continue

                                    # Check for internal introns
                                    firstBreakFound = False
                                    secondBreakFound = False
                                    # Check in data
                                    for this_start in all_coord_start2end[chr]:
                                        if exon_start < this_start < exon_end:
                                            firstBreakFound = True
                                            # Reduce search time if second
                                            # break is same intron
                                            for this_end in all_coord_start2end[chr][this_start]:
                                                if this_end < exon_end:
                                                    secondBreakFound = True
                                                    break
                                            break

                                    if firstBreakFound and not secondBreakFound:
                                        for this_end in all_coord_end2start[chr]: 
                                            if exon_start < this_end < exon_end:
                                                # This has to be a full internal
                                                # intron
                                                if this_start < this_end:
                                                    secondBreakFound = True
                                                    break

                                    if firstBreakFound and secondBreakFound:
                                        continue

                                    # Check in annotation
#                                    if hasInternalIntron(coord_start2end,
#                                                         chr,
#                                                         exon_start, exon_end):
#                                        continue

                                    if hasInternalTerminationInitationExon(alt_first_exons,
                                                                           chr,
                                                                           exon_start,
                                                                           exon_end):
                                        continue
                                    if hasInternalTerminationInitationExon(alt_last_exons,
                                                                           chr,
                                                                           exon_start,
                                                                           exon_end):
                                        continue
                                        
                                    # Add to excl junction set
                                    if chr in excl_jcns:
                                        excl_jcns[chr].add((start, end))
                                    else:
                                        excl_jcns[chr] = set([(start, end)])                                      
 
                                    # Found a cassette event        
                                    if exon_coord in cassette_exon_dict:
                                        cassette_exon_dict[exon_coord]["excl_set"].add(excl_str)
                                        cassette_exon_dict[exon_coord]["left_set"].add(left_str)
                                        cassette_exon_dict[exon_coord]["right_set"].add(right_str)
                                        cassette_exon_dict[exon_coord]["strand"] = strand
                                    else: 
                                        cassette_exon_dict[exon_coord] = {"excl_set": set([excl_str]),
                                                                        "left_set": set([left_str]),
                                                                        "right_set": set([right_str]),
                                                                        "strand": strand}

    for exon_coord in cassette_exon_dict:

        chr = exon_coord[0]
        exon_start = exon_coord[1]
        exon_end = exon_coord[2]

        isNovel = False

        # Print novel cassette
        excl_raw_count = 0
        incl_raw_count = 0
        excl_lenNorm_count = 0
        incl_lenNorm_count = 0

        if full_exon_count_dict:
            # Inclusion counts only come from exon counts
            exon_str = formatCoordStr(chr, exon_start, exon_end)
            exonFound = False
#           if exon_str in full_exon_count_dict:
#               exonFound = True
#               incl_file1_count += full_exon_count_dict[exon_str][0] 
#               incl_file2_count += full_exon_count_dict[exon_str][1] 
#           if exon_str in full_multi_exon_count_dict:
#               exonFound = True
#               incl_file1_count += full_multi_exon_count_dict[exon_str][0] 
#               incl_file2_count += full_multi_exon_count_dict[exon_str][1] 
#           if not exonFound:
#               print "No counts for cassette exon: %s" % exon_str
#           
#           # Paired end counting for constitutive portion
# BOOKMARK!!!!!
#            if 
#                       start_multi_exon_count_dict,
#                       end_multi_exon_count_dict,
#            def findAdjacentSharedRegion(chr, annotated_exons, pos, n_or_p):

        else:
            left_jcn_strs = []
            left_jcn_counts_raw = []
#            left_jcn_counts_lenNorm = []

            right_jcn_strs = []
            right_jcn_counts_raw = []
#            right_jcn_counts_lenNorm = []

            for left_str in cassette_exon_dict[exon_coord]["left_set"]:
                left_str_raw_ct = all_jcn_count_dict[left_str][1]

                if norm2:
                    left_str_raw_ct = int(round(left_str_raw_ct/norm2))

                incl_raw_count += left_str_raw_ct

                if not isAnnotated(left_str, annotated_introns):
                    isNovel = True

                left_jcn_strs.append(left_str)
                left_jcn_counts_raw.append(left_str_raw_ct)

            for right_str in cassette_exon_dict[exon_coord]["right_set"]:
                right_str_raw_ct = all_jcn_count_dict[right_str][1]

                if norm2:
                    right_str_raw_ct = int(round(right_str_raw_ct/norm2))

                incl_raw_count += right_str_raw_ct

                if not isAnnotated(right_str, annotated_introns):
                    isNovel = True

                right_jcn_strs.append(right_str)
                right_jcn_counts_raw.append(right_str_raw_ct)

        excl_jcn_strs = []
        excl_jcn_counts_raw = []
        # For both paired-end counting and single end counting, the exclusion
        # counts come from the exclusion junction
        left_most_excl_start = INFINITY
        right_most_excl_end = 0
        for excl_str in cassette_exon_dict[exon_coord]["excl_set"]:
            excl_str_raw_ct = all_jcn_count_dict[excl_str][1]

            if norm2:
                excl_str_raw_ct = int(round(excl_str_raw_ct/norm2))

            excl_raw_count += excl_str_raw_ct
        
            if not isAnnotated(excl_str, annotated_introns):
                isNovel = True

            excl_chr, excl_start, excl_end = convertCoordStr(excl_str)
            if excl_start < left_most_excl_start:
                left_most_excl_start = excl_start
            if excl_end > right_most_excl_end:
                right_most_excl_end = excl_end

            excl_jcn_strs.append(excl_str)
            excl_jcn_counts_raw.append(excl_str_raw_ct)

        # Find flanking constitutive regions
        upstr_const = findAdjacentSharedRegion(chr,
                                               cassette_exon_dict[exon_coord]["strand"],
                                               annotated_exons_by_strand, left_most_excl_start - 1, "P")

        dwnstr_const = findAdjacentSharedRegion(chr, 
                                                cassette_exon_dict[exon_coord]["strand"],
                                                annotated_exons_by_strand, right_most_excl_end + 1, "N")

        const_strs = []

        if upstr_const:
            const_strs.append(upstr_const)
        if dwnstr_const:
            const_strs.append(dwnstr_const)
            

        # Add these regions to exon_coords for quantification:
        if printExonCoord:
            if upstr_const:
                exon_coords.add(convertCoordStr(upstr_const))
            if dwnstr_const:
                exon_coords.add(convertCoordStr(dwnstr_const))

        # A cassette exon is unannotated if any junction is novel.
        label = "K"
        if isNovel:
            label = "N"

        # Normalize counts by length. Even though multiple junctions are
        # involved in a given inclusion or exclusion event, the normalization
        # factor is only on one junction
        exclusion_len = jcn_seq_len
        inclusion_len = (2* jcn_seq_len) + exon_end - exon_start + 1

        excl_lenNorm_count = normalizeByLen(excl_raw_count, exclusion_len)
        incl_lenNorm_count = normalizeByLen(incl_raw_count, inclusion_len)

        # Add Exclusion or Inclusion Annotation
        # Checks percent exclusion from both files.
#       e_or_i = checkExclusionInclusion(excl_file1_count,
#                                       incl_file1_count,
#                                       excl_file2_count,
#                                       incl_file2_count) 

        gene_name = inferGeneName(annotated_genes_by_strand, chr, cassette_exon_dict[exon_coord]["strand"],
                                  [(exon_start, exon_end)])


        # print
        out_str = "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d" % (label,
                        "?",
                        gene_name,
                        chr,
                        cassette_exon_dict[exon_coord]["strand"],
                        ",".join(list(cassette_exon_dict[exon_coord]["left_set"])),
                        exon_start,
                        exon_end,
                        ",".join(list(cassette_exon_dict[exon_coord]["right_set"])),
                        excl_raw_count,
                        incl_raw_count,
                        excl_lenNorm_count,
                        incl_lenNorm_count)                                                

        if printExonCoord:
            ce_str = "\t%s" % formatCoordStr(chr, exon_start, exon_end)
            out_str += ce_str
            
            exon_coords.add((chr, exon_start, exon_end))

        cassette_out.write(out_str + "\n")
        if chr in return_cassette_exons:
            return_cassette_exons[chr].add((exon_start, exon_end))
        else:
            return_cassette_exons[chr] = set([(exon_start, exon_end)])

        # Print to all_event_file
        out_str = getAllEventStr(label,
                              "cassette",
                              gene_name,
                              chr,
                              cassette_exon_dict[exon_coord]["strand"],
                              ";".join(excl_jcn_strs),
                              ",".join(left_jcn_strs)+ ";" + ",".join(right_jcn_strs), 
                              "",
                              "%s" % formatCoordStr(chr, exon_start, exon_end), 
                              "",
                              ";".join(const_strs),
                              ";".join(map(repr,excl_jcn_counts_raw)), 
                              ",".join(map(repr,left_jcn_counts_raw)) + ";" + ",".join(map(repr,right_jcn_counts_raw)),
                              repr(excl_raw_count),
                              repr(incl_raw_count),
                              "", "", 
                              None, None, 
                              "",
                              None, 
                              "",
                              None)
       
        all_event_info_out.write(out_str + "\n")
        
    return return_cassette_exons 

def printIREvents(db, annotated_genes, annotated_genes_by_strand, annotated_exons,
                  annotated_exons_by_strand, coord_start2end, coord_end2start,
                  jcn_count_dict, jcn2strand, ir_count_dict, 
                  ir_left_out, ir_right_out, printExonCoords, exon_coords,
                  norm1, norm2, jcn_seq_len):
    """
    Print out all events
    """
    # Print left
    for chr in coord_start2end:
        for start in coord_start2end[chr]:

            excl_file_raw_count = 0
            incl_file_raw_count = 0
            excl_file_lenNorm_count = 0
            incl_file_lenNorm_count = 0
        
            jcn_str_list = []

            strand = None

            good_jcn_str = None

            for end in coord_start2end[chr][start]:

                jcn_str = formatCoordStr(chr, start, end)
                jcn_str_list.append(jcn_str)
        
                strand = updateStrand(strand, jcn2strand[jcn_str])

                if jcn_str not in jcn_count_dict or \
                   jcn_str not in ir_count_dict:
                    continue
        
                # jcn_str exists in ir_count_dict
                good_jcn_str = jcn_str

#               excl_file1_count += normalizeByLen(jcn_count_dict[jcn_str][0],
#                                                  jcn_seq_len)
                excl_file_raw_count += jcn_count_dict[jcn_str][1]
                excl_file_lenNorm_count += normalizeByLen(jcn_count_dict[jcn_str][1],
                                                          jcn_seq_len)

            if good_jcn_str:
                # Only need to add inclusion counts from one of the introns
#                incl_file1_count += ir_count_dict[good_jcn_str]["left"][0]
                incl_file_raw_count += ir_count_dict[good_jcn_str]["left"][1]

#            if excl_file1_count == 0 and excl_file2_count == 0 and incl_file1_count == 0 and incl_file2_count == 0:
#                continue

            # Add potential constitutive exon to exon_coords
            upstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, start - 1, "P")
            if printExonCoords:
                if upstr_const:
                    exon_coords.add(convertCoordStr(upstr_const))
                

            # Check if it is a novel junction or not
#                label = "K"
#            else:
#                label = "N"

            gene_name = inferGeneName(annotated_genes_by_strand, chr, strand,
                                      [(start, end)])

            if norm2:
#               excl_file1_count = int(round(excl_file1_count/norm1))
#               incl_file1_count = int(round(incl_file1_count/norm1))
            
                excl_file_raw_count = int(round(excl_file_raw_count/norm2))
                incl_file_raw_count = int(round(incl_file_raw_count/norm2))

#            incl_file1_count = normalizeByLen(incl_file1_count, jcn_seq_len)
            incl_file_lenNorm_count = normalizeByLen(incl_file_raw_count, jcn_seq_len)

            # Checks percent exclusion from both files.
#           e_or_i = checkExclusionInclusion(excl_file1_count,
#                                           incl_file1_count,
#                                           excl_file2_count,
#                                           incl_file2_count) 

            out_str = "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n" % ("?",
                                                            gene_name,
                                                            chr,
                                                            strand,
                                                            start,
                                                            ",".join(jcn_str_list),
                                                            excl_file_raw_count,
                                                            incl_file_raw_count,
                                                            excl_file_lenNorm_count,
                                                            incl_file_lenNorm_count)

            if hasNegativeVals(incl_file_raw_count, incl_file_lenNorm_count, 0, 0):
                ERROR_LOG.write("Negative Vals: %s\n" % out_str)
                out_str = "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n" % ("?",
                                                            gene_name,
                                                            chr,
                                                            strand,
                                                            start,
                                                            ",".join(jcn_str_list),
                                                            0,0,0,0)

            ir_left_out.write(out_str)

    # Print right
    for chr in coord_end2start:
        for end in coord_end2start[chr]:

            excl_file_raw_count = 0
            incl_file_raw_count = 0
            excl_file_lenNorm_count = 0
            incl_file_lenNorm_count = 0
        
            jcn_str_list = []

            strand = None

            good_jcn_str = None

            for start in coord_end2start[chr][end]:

                jcn_str = formatCoordStr(chr, start, end)
                jcn_str_list.append(jcn_str)

                strand = updateStrand(strand, jcn2strand[jcn_str])

                if jcn_str not in jcn_count_dict or \
                   jcn_str not in ir_count_dict:
                    continue

                # Good junction string will exist in the ir dictionary
                good_jcn_str = jcn_str

#                excl_file1_count += jcn_count_dict[jcn_str][0]
                excl_file_raw_count += jcn_count_dict[jcn_str][1]

            if good_jcn_str:
                # Only need to add inclusion counts from one of the introns
#                incl_file1_count += ir_count_dict[good_jcn_str]["right"][0]
                incl_file_raw_count += ir_count_dict[good_jcn_str]["right"][1]
                    
#            if excl_file1_count == 0 and excl_file2_count == 0 and incl_file1_count == 0 and incl_file2_count == 0:
#                continue

            downstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, end + 1, "N")
            if printExonCoords:
                if downstr_const:
                    exon_coords.add(convertCoordStr(downstr_const))
 

            # Check if it is a novel ir event or not
#            if jcn_str in ir_introns:
#                label = "K"
#            else:
#                label = "N"

            
            gene_name = inferGeneName(annotated_genes_by_strand, chr, strand,
                                      [(start,end)])
    
            if norm2:
#               excl_file1_count = int(round(excl_file1_count/norm1))
#               incl_file1_count = int(round(incl_file1_count/norm1))
            
                excl_file_raw_count = int(round(excl_file_raw_count/norm2))
                incl_file_raw_count = int(round(incl_file_raw_count/norm2))

#           excl_file1_count = normalizeByLen(excl_file1_count, 
#                                             jcn_seq_len)
            excl_file_lenNorm_count = normalizeByLen(excl_file_raw_count,
                                              jcn_seq_len)
            
#            incl_file1_count = normalizeByLen(incl_file1_count, jcn_seq_len)
            incl_file_lenNorm_count = normalizeByLen(incl_file_raw_count, jcn_seq_len)

            # Checks percent exclusion from both files.
#           e_or_i = checkExclusionInclusion(excl_file1_count,
#                                           incl_file1_count,
#                                           excl_file2_count,
#                                           incl_file2_count) 

            out_str = "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n" % ("?",
                                                            gene_name,
                                                            chr,
                                                            strand,
                                                            end,
                                                            ",".join(jcn_str_list),
                                                            excl_file_raw_count,
                                                            incl_file_raw_count,
                                                            excl_file_lenNorm_count,
                                                            incl_file_lenNorm_count)

            if hasNegativeVals(incl_file_raw_count, incl_file_lenNorm_count, 0, 0):
                ERROR_LOG.write("Negative Vals: %s\n" % out_str)
                out_str = "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n" % ("?",
                                                            gene_name,
                                                            chr,
                                                            strand,
                                                            start,
                                                            ",".join(jcn_str_list),
                                                            0,0,0,0)

            ir_right_out.write(out_str)

def printMultiCassetteExons(db,
                            annotated_genes,
                            annotated_genes_by_strand,
                            alt_first_exons,
                            alt_last_exons,
                            alt_first_exon_search_tree,
                            alt_last_exon_search_tree,
                            all_jcn_count_dict, 
                            all_coord_start2end,
                            all_coord_end2start,
                            all_jcn2strand,
                            all_jcn_search_tree,
                            full_exon_count_dict,
                            full_multi_exon_count_dict,
                            cassette_exons,
                            annotated_introns,
                            annotated_exons,
                            annotated_exons_by_strand,
                            printExonCoords,
                            exon_coords,
                            excl_jcns,
                            mc_out,
                            all_event_info_out,
                            norm1, norm2, jcn_seq_len):
    """
    Some of the novel cassette exons will pick up multi-cassette events.  I am
    listing those as well.

    The alt_first_exons and alt_last_exons dictionary is for some of the multi
    cassette events that are miscategorized because of these events.
    """
    for chr in all_coord_start2end:
        for exclusion_start in all_coord_start2end[chr]:
            if len(all_coord_start2end[chr][exclusion_start]) < 2:
                continue
            for exclusion_end in all_coord_start2end[chr][exclusion_start]:
                if len(all_coord_end2start[chr][exclusion_end]) < 2:
                    continue
                for left_inclusion_end in all_coord_start2end[chr][exclusion_start]:
                    if left_inclusion_end == exclusion_end:
                        continue
                    for right_inclusion_start in all_coord_end2start[chr][exclusion_end]:
                        if right_inclusion_start <= left_inclusion_end:
                            continue
                        # Check for potential internal cassette exons by
                        # looking for internal introns that exist between the
                        # left inclusion end and the right inclusion_start
                        excl_jcn = formatCoordStr(chr, 
                                                 exclusion_start,
                                                 exclusion_end)
                        strand = None
                        strand = updateStrand(strand, all_jcn2strand[excl_jcn])
        
                        internal_introns = getInternalIntrons(all_jcn_search_tree[chr][strand],
                                                              left_inclusion_end,
                                                              right_inclusion_start)

                        if internal_introns is None:
                            continue

                        internal_introns_len = len(internal_introns)
                        if internal_introns_len > 1:
                            if internal_introns_len > MAX_EXON_CLUSTER:
#                                print "Multi-cassette: Exon cluster is too large=%d" % internal_introns_len
                                continue
                            intron_clusters = findNonOverlappingSets(-1,
                                                                     internal_introns)
                            intron_clusters = removeRedundantSets(intron_clusters)
                        else:
                            intron_clusters = [set(internal_introns)]
        
                        if internal_introns == []:
                            continue

                        # Each cluster is a potential multi-cassette event
                        for intron_cluster in intron_clusters:
                            intron_list = list(intron_cluster)
                            intron_list.sort()

                            # Check for novel introns
                            isNovel = False

                            if chr not in annotated_introns:
                                isNovel = True
                            elif ((exclusion_start, exclusion_end, "+") not in annotated_introns[chr]) and \
                               ((exclusion_start, exclusion_end, "-") not in annotated_introns[chr]):
                                isNovel = True
                            elif ((exclusion_start, left_inclusion_end, "+") not in annotated_introns[chr]) and \
                                 ((exclusion_start, left_inclusion_end, "-") not in annotated_introns[chr]):
                                isNovel = True
                            elif ((right_inclusion_start, exclusion_end, "+") not in annotated_introns[chr]) and \
                               ((right_inclusion_start, exclusion_end, "-") not in annotated_introns[chr]):
                                isNovel = True
                        

                            mc_exon_list = []
                            prev_intron_end = left_inclusion_end
                            for intron_coord in intron_list:
                                mc_exon_list.append((prev_intron_end + 1,
                                                     intron_coord[0] - 1))
                                prev_intron_end = intron_coord[1] 

                                # Check for novel introns
                                if chr not in annotated_introns:
                                    isNovel = True
                                elif ((intron_coord[0], intron_coord[1], "+") not in annotated_introns[chr]) and \
                                   ((intron_coord[0], intron_coord[1], "-") not in annotated_introns[chr]):
                                    isNovel=True

                            # Add last exon
                            mc_exon_list.append((prev_intron_end + 1,
                                                 right_inclusion_start - 1))

                            # Double check that all exons in the list have no
                            # internal introns and are not altnernative first
                            # and last exons
                            exonCheck = True
                            for (exon_start, exon_end) in mc_exon_list:
                                if hasInternalIntron(all_coord_start2end, chr,
                                                     exon_start, exon_end):
                                    exonCheck = False
                                    break
            
                                a_exonFound = False
                                if annotated_exons:
                                    if chr in annotated_exons:
                                        if (exon_start, exon_end, "+") in annotated_exons[chr]:
                                            a_exonFound = True
                                        if (exon_start, exon_end, "-") in annotated_exons[chr]:
                                            a_exonFound = True

                                    if not a_exonFound:
                                        exonCheck = False
                                        break

                                if hasOverlap(alt_first_exon_search_tree[chr],
                                              (exon_start, exon_end)):
                                    exonCheck=False
                                    break

                                if hasOverlap(alt_last_exon_search_tree[chr],
                                              (exon_start, exon_end)):
                                    exonCheck=False
                                    break
           
                                if not exonCheck:
                                    break

                            if not exonCheck:
                                continue

                            # Format exons for output
                            exon_strs = []
                            len_of_exons = 0
                            for (exon_start, exon_end) in mc_exon_list:
                                # Add to exon_strs
                                exon_str = formatCoordStr(chr, exon_start, exon_end)
                                exon_strs.append(exon_str)
                                len_of_exons += (exon_end - exon_start + 1)
                                # Add to cassette exon dict
                                if chr in cassette_exons:
                                    cassette_exons[chr].add((exon_start, exon_end))
                                else:
                                    cassette_exons[chr] = set([(exon_start,
                                                                exon_end)])
                            # Multi-cassette event found
                            excl_file_raw_count = 0
                            incl_file_raw_count = 0
                            excl_file_lenNorm_count = 0
                            incl_file_lenNorm_count = 0
    
                            # Add to exclusion junctions
                            if chr in excl_jcns:
                                excl_jcns[chr].add((exclusion_start,    
                                                    exclusion_end))
                            else:   
                                excl_jcns[chr] = set([(exclusion_start,
                                                       exclusion_end)])

                            # Exclusion_count
                            excl_jcn_str = formatCoordStr(chr, exclusion_start,
                                                         exclusion_end)
                            excl_file_raw_count = all_jcn_count_dict[excl_jcn_str][1]

                            if norm2:
                                excl_file_raw_count = int(round(all_jcn_count_dict[excl_jcn_str][1]/norm2))

                            strand = None
                            strand = updateStrand(strand, all_jcn2strand[excl_jcn_str])

                            upstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, exclusion_start - 1, "P")
                            dwnstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, exclusion_end + 1, "N")

                            const_strs = []
                            if upstr_const:
                                const_strs.append(upstr_const)
                                if printExonCoords:
                                    exon_coords.add(convertCoordStr(upstr_const))
                            if dwnstr_const:
                                const_strs.append(dwnstr_const)
                                if printExonCoords:
                                    exon_coords.add(convertCoordStr(dwnstr_const))
                            
                            inclusion_jcns = [] 
                            raw_inclusion_cts = [] # parallel to inclusion_jcns [raw_ct]
                            lenNorm_inclusion_cts = [] # parallel to inclusion_jcns [lenNorm_ct]

                            # Paired end counting needs to be added
                            if full_exon_count_dict:
                                for exon_str in exon_strs:
                                    exonFound = False
                                    if exon_str in full_exon_count_dict:
                                        exonFound = True
                                        incl_file1_count += \
                                                full_exon_count_dict[exon_str][0]
                                        incl_file2_count += \
                                                full_exon_count_dict[exon_str][1]
                                    if exon_str in full_multi_exon_count_dict:
                                        exonFound = True
                                        incl_file1_count += \
                                                full_multi_exon_count_dict[exon_str][0]
                                        incl_file2_count += \
                                                full_multi_exon_count_dict[exon_str][1]
                                    if not exonFound:
                                        print "No counts for multi-cassette: %s" %\
                                               exon_str
                        
                            else: # Just single reads, add junction counts
                                # Inclusion_count
                                last_start = exclusion_start
                                last_end = left_inclusion_end
                                last_jcn = formatCoordStr(chr, last_start, last_end)

                                incl_file_raw_count = 0


                                for intron_coord in intron_list:
                                    this_start = intron_coord[0]
                                    this_end = intron_coord[1]
                                    this_jcn = formatCoordStr(chr, this_start, this_end)

                                    inferred_exon = formatCoordStr(chr,
                                                                  last_end + 1,
                                                                  this_start - 1)

                                    (mc_proportion1,
                                     mc_proportion2) = getMXE_MCProportion(inferred_exon,
                                                                           last_start,
                                                                           this_end,
                                                                           all_jcn_count_dict,
                                                                           all_coord_start2end,
                                                                           all_coord_end2start)
                
                                    # proportions used to estimate junction count for
                                    # mutually exclusive events based on all junctions on
                                    # each side of the exon
                                    upstrm_jcn_sum_raw = 0

#                                   The end coordinate shouldn't have been
#                                   this_end bu last end
#                                    for other_start in all_coord_end2start[chr][this_end]:
                                    for other_start in all_coord_end2start[chr][last_end]:
#                                       upstrm_jcn_sum1 += all_jcn_count_dict["%s_%d_%d" % (chr,
#                                                                                          other_start,
#                                                                                          this_end)][0]
#                                       upstrm_jcn_sum_raw += all_jcn_count_dict[formatCoordStr(chr,
#                                                                                          other_start,
#                                                                                          this_end)][1]
                                        upstrm_jcn_sum_raw += all_jcn_count_dict[formatCoordStr(chr,
                                                                                           other_start,
                                                                                           last_end)][1]
                            
#                                   this_incl_file1_count = int(round(upstrm_jcn_sum1 *
#                                                                     mc_proportion1))
                                                            
                                    this_incl_file_raw_count = int(round(upstrm_jcn_sum_raw *
                                                                      mc_proportion2))

                                    if norm2:
                                        this_incl_file_raw_count = int(round(this_incl_file_raw_count/norm2))

#                                    incl_file1_count += this_incl_file1_count
                                    incl_file_raw_count += this_incl_file_raw_count

                                    inclusion_jcns.append(last_jcn)
                                    raw_inclusion_cts.append(this_incl_file_raw_count)

                                    # Now update variables
                                    last_start = this_start
                                    last_end = this_end
                                    last_jcn = this_jcn


                                # Now add junctions from last exon
                                inferred_exon = formatCoordStr(chr, last_end + 1,
                                                              right_inclusion_start - 1)
                                                           
                                (mc_proportion1,
                                 mc_proportion2) = getMXE_MCProportion(inferred_exon, 
                                                                       last_start,
                                                                       exclusion_end,
                                                                       all_jcn_count_dict,
                                                                       all_coord_start2end,
                                                                       all_coord_end2start)

                                upstrm_jcn = formatCoordStr(chr, last_start,
                                                           last_end)
                                dwnstrm_jcn = formatCoordStr(chr, right_inclusion_start,
                                                            exclusion_end)

                                # proportions used to estimate junction count for
                                # mutually exclusive events based on all junctions on
                                # each side of the exon
                                upstrm_jcn_sum_raw = 0
                                dwnstrm_jcn_sum_raw = 0

                                for other_start in all_coord_end2start[chr][last_end]:
#                                   upstrm_jcn_sum1 += all_jcn_count_dict["%s_%d_%d" % (chr,
#                                                                                      other_start,
#                                                                                      last_end)][0]
                                    upstrm_jcn_sum_raw += all_jcn_count_dict[formatCoordStr(chr,
                                                                                       other_start,
                                                                                       last_end)][1]
                                for other_end in all_coord_start2end[chr][right_inclusion_start]:
#                                   dwnstrm_jcn_sum1 += all_jcn_count_dict["%s_%d_%d" % (chr,
#                                                                                        right_inclusion_start,
#                                                                                        other_end)][0]
                                    dwnstrm_jcn_sum_raw += all_jcn_count_dict[formatCoordStr(chr,
                                                                                         right_inclusion_start,
                                                                                         other_end)][1]

#                               upstrm_incl_file1_count = int(round(upstrm_jcn_sum1 *
#                                                                   mc_proportion1))
                                upstrm_incl_file_raw_count = int(round(upstrm_jcn_sum_raw *
                                                                    mc_proportion2))
#                               dwnstrm_incl_file1_count = int(round(dwnstrm_jcn_sum1 *
#                                                                   mc_proportion1))
                                dwnstrm_incl_file_raw_count = int(round(dwnstrm_jcn_sum_raw*
                                                                    mc_proportion2))

                                if norm2:
                                    upstrm_incl_file_raw_count = int(round(upstrm_incl_file_raw_count/norm2))
                                    dwnstrm_incl_file_raw_count = int(round(dwnstrm_incl_file_raw_count/norm2))

#                                incl_file1_count += upstrm_incl_file1_count
                                incl_file_raw_count += upstrm_incl_file_raw_count
#                                incl_file1_count += dwnstrm_incl_file1_count
                                incl_file_raw_count += dwnstrm_incl_file_raw_count
                
                                inclusion_jcns.append(upstrm_jcn)
                                inclusion_jcns.append(dwnstrm_jcn)
                                raw_inclusion_cts.append(upstrm_incl_file_raw_count)
                                raw_inclusion_cts.append(dwnstrm_incl_file_raw_count)

                            # Calculate inclusion length for normalization       
                            incl_length = (len(inclusion_jcns)*jcn_seq_len) + len_of_exons 
 
                            excl_file_lenNorm_count = normalizeByLen(excl_file_raw_count, jcn_seq_len)

                            incl_file_lenNorm_count = normalizeByLen(incl_file_raw_count, incl_length)

                            # Add Exclusion or Inclusion Annotation
                            # Checks percent exclusion from both files.
#                           e_or_i = checkExclusionInclusion(excl_file1_count,
#                                                           incl_file1_count,
#                                                           excl_file2_count,
#                                                           incl_file2_count) 

                            label = "K"
                            if isNovel:
                                label = "N"

                            gene_name = inferGeneName(annotated_genes_by_strand, chr, strand,
                                                      inclusion_jcns + [excl_jcn_str])


                            # print
                            out_str = "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d" % (label,
                                                                                    "?",
                                                                                    gene_name,
                                                                                    chr,
                                                                                    strand,
                                                                                    exclusion_start,
                                                                                    exclusion_end,
                                                                                    ",".join(exon_strs),
                                                                                    excl_file_raw_count,
                                                                                    incl_file_raw_count,
                                                                                    excl_file_lenNorm_count,
                                                                                    incl_file_lenNorm_count)                                                
                            if printExonCoords:
                                out_str += "\t%s" % ",".join(exon_strs)
            
                                for exon_s in exon_strs:
                                    exon_coords.add(convertCoordStr(exon_s))

                            mc_out.write(out_str + "\n")

                            # Now print to all event file
                            out_str = getAllEventStr(label, "coord_cassette", gene_name,  chr, strand, 
                                                     excl_jcn_str, ";".join(inclusion_jcns),
                                                     "", ";".join(exon_strs),
                                                     "",
                                                     ";".join(const_strs),
                                                     excl_file_raw_count, ";".join(map(repr,raw_inclusion_cts)),
                                                     repr(excl_file_raw_count),
                                                     repr(incl_file_raw_count),
                                                     "","",
                                                     None,None,
                                                     "",
                                                     None,
                                                     "",
                                                     None)                                                

                            all_event_info_out.write(out_str + "\n")


def printMutuallyExclusive(db,
                           annotated_genes,
                           annotated_genes_by_strand,
                           annotated_internal_exons,
                           alt_first_exons,
                           alt_last_exons,
                           all_jcn_count_dict,
                           all_coord_start2end,
                           all_coord_end2start,
                           all_jcn2strand,
                           full_exon_count_dict,
                           full_multi_exon_count_dict,
                           cassette_exons,
                           annotated_introns,
                           annotated_exons,
                           annotated_exons_by_strand,
                           annotated_exon_search_tree,
                           printExonCoords,
                           exon_coords,
                           excl_jcns,
                           me_out,
                           all_event_info_out,
                           norm1, norm2,
                           jcn_seq_len):
    """
    Mutually exclusive exons are a set of cassette exons whose upstream
    introns start at the same position and the downstream introns end at the
    same position.  There should be no internal introns within the exonic
    regions.
    """
    for chr in all_coord_start2end:
        print "Mutually exclusive %s" % chr
        # Create mutually exclusive dictionary to look through
        # The dictionary is of the form:
        # {(upstream intron start, downstream intron end):[[(exon_start, exon_end),], 
        #                                                  [(exon_start, exon_end),],]}
        me_dict2 = buildMutuallyExclusiveDict(chr, all_coord_start2end, all_coord_end2start)
            
        # me_dict2 is a result of filtering me_dict
#       me_dict2 = {}
#       # Check for internal introns in the dictionary and remove those clusters.
#       for me_coord in me_dict.keys():
#           for exon_list in me_dict[me_coord]:    
#               new_exon_list = []
#               coord_copy = list(exon_list)
#               for (exon_start, exon_end) in coord_copy:
#                   if not hasInternalIntron(all_coord_start2end, chr, exon_start,
#                                        exon_end):
#                       new_exon_list.append((exon_start, exon_end)) 
#                        me_dict[me_coord][exon_list_idx].remove((exon_start, exon_end))
# Not sure why this break statement is here...
#                    break # to next me_coord

                # if all but one or no exons are left, then remove the
                # exon_list
#               if len(new_exon_list) >= 2:
#                   updateDictOfLists(me_dict2, me_coord, new_exon_list)

        # Check that all are annotated exons
        if annotated_exons:
            me_dict3 = {}
            for me_coord in me_dict2.keys():
                for exon_list in me_dict2[me_coord]:
                    new_exon_list = []
                    coord_copy = list(exon_list)
                    for (exon_start, exon_end) in coord_copy:
                        exonFound = False
                        if chr in annotated_exons:
                            if (exon_start, exon_end, "+") in annotated_exons[chr]:
                                exonFound = True
                            elif (exon_start, exon_end, "-") in annotated_exons[chr]:
                                exonFound = True
                        if exonFound:
                            new_exon_list.append((exon_start, exon_end))
#                            me_dict[me_coord].remove((exon_start, exon_end))

                    if len(new_exon_list) >= 2:
                        updateDictOfLists(me_dict3, me_coord, new_exon_list)

        else:
            me_dict3 = dict(me_dict2)
            
        # NO LONGER ENFORCED 110223: Make sure that no connection connects any of the exons in the group
#       me_dict4 = {}
#       for me_coord in me_dict3.keys():
#           for exon_list in me_dict3[me_coord]:
#               if not isConnected(exon_list, chr, all_coord_start2end):
#                   updateDictOfSets(me_dict4, me_coord, tuple(exon_list))

        # Make sure that there are no other annotated exons intervening the
        # mutually exclusive exon cluster
        if annotated_exons:
            me_dict4 = {}

            for me_coord in me_dict3:           
                upstream_start = me_coord[0]
                downstream_end = me_coord[1]
                for exon_list in me_dict3[me_coord]:
                    list_copy = list(exon_list)
                    # Sort the list of exons
                    list_copy.sort()
                    intervening_spaces = []  
                    strand = None
                    last_end = list_copy[0][1]
                    for (exon_start, exon_end) in list_copy[1:]:
                        # Infer strand information
                        jcn = formatCoordStr(chr, upstream_start, exon_start - 1)
                        strand = updateStrand(strand, all_jcn2strand[jcn])
                        intervening_spaces.append((last_end + 1, exon_start - 1))
                        last_end = exon_end

                    # Now check for any completely contained exons in any of the
                    # intervening spaces
                    internalExonFound = False
                    if chr in annotated_exon_search_tree:
                        if strand in annotated_exon_search_tree[chr]:
                            for (space_start, space_end) in intervening_spaces:
                                if hasOuterContainer(annotated_exon_search_tree[chr][strand],
                                                     (space_start, space_end)):
                                    internalExonFound = True
                                    break
#                           for (exon_start, exon_end, exon_strand) in annotated_exons[chr]:
#                               if exon_strand != strand:
#                                   continue

#                               if (exon_start, exon_end) in alt_first_exons[chr]:
#                                   continue
#                               if (exon_start, exon_end) in alt_last_exons[chr]:
#                                   continue
#                               if (exon_start, exon_end) not in annotated_internal_exons[chr]:
#                                   continue

#                               if space_start < exon_start < space_end:
#                                   if space_start < exon_end < space_end:
#                                       internalExonFound = True
#                                       break
#                           if internalExonFound:
#                               break

                    if not internalExonFound:
                        updateDictOfLists(me_dict4, me_coord, list_copy)
                            
        else:
            me_dict4 = dict(me_dict3)

        # Now, calculate inclusion and exclusion for each exon in ME events
        for (upstream_start, downstream_end) in me_dict4:
            # Infer strand information
            strand = None
            jcn = formatCoordStr(chr, upstream_start, 
                                me_dict4[(upstream_start, downstream_end)][0][0][0] - 1)  # one of the exon's start -1

            strand = updateStrand(strand, all_jcn2strand[jcn])

            upstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, upstream_start - 1, "P")
            dwnstr_const = findAdjacentSharedRegion(chr, strand, annotated_exons_by_strand, downstream_end + 1, "N")

            const_strs = []

            if upstr_const:
                const_strs.append(upstr_const)
            if dwnstr_const:
                const_strs.append(dwnstr_const)
                
            # Add these regions to exon_coords for quantification:
            if printExonCoords:
                if upstr_const:
                    exon_coords.add(convertCoordStr(upstr_const))
                if dwnstr_const:
                    exon_coords.add(convertCoordStr(dwnstr_const))

            for exon_list in me_dict4[(upstream_start, downstream_end)]:
                total_me_cts = [0,0]
                inclusion_cts = [] # parallel to me_dict exons
                exon_strs = [] # parallel to me_dict exons for the strings of the
                               # exons
                detailed_jcns = [] # Parallel to exon strs
                                   # [(left_jcn, right_jcn),]
                detailed_inclusion_cts = [] # parallel to inclusion_cts but
                                            # gives information for both
                                            # junctions
                                            # [((left_jcn_raw, left_jcn_lenNorm),(right_jcn_raw, right_jcn_lenNorm)),]
                isNovel = False # If any of the introns are novel = True
                # Get inclusion counts for each exon
        
                strand = None

                for (exon_start, exon_end) in exon_list:

                    # Add to exon_strs
                    exon_str = formatCoordStr(chr, exon_start, exon_end)
                    exon_strs.append(exon_str)


                    # Update strand information for both upstream and
                    # downstream introns
                    upstr_str = formatCoordStr(chr, upstream_start, exon_start - 1)
                    dwnstr_str = formatCoordStr(chr, exon_end + 1, downstream_end)

                    strand = updateStrand(strand, all_jcn2strand[upstr_str])
                    strand = updateStrand(strand, all_jcn2strand[dwnstr_str])

                    # Count paired end
                    if full_exon_count_dict:
                        # Add exon counts to total
                        exonFound = False
#                       if full_exon_count_dict is not None:
#                           if exon_str in full_exon_count_dict:
#                               exonFound = True
#                               raw_ct = full_exon_count_dict[exon_str][0]
#                               lenNorm_ct = full_exon_count_dict[exon_str][1]
#                           if exon_str in full_multi_exon_count_dict:
#                               exonFound = True
#                               raw_ct += full_multi_exon_count_dict[exon_str][0]
#                               lenNorm_ct += full_multi_exon_count_dict[exon_str][1]

#                           if not exonFound:
#                               print "No counts for mutually exclusive exon: %s" %\
#                                      exon_str

                    else: # Just add junction counts

                        left_jcn_str = formatCoordStr(chr, upstream_start, 
                                                     exon_start - 1)
                        right_jcn_str = formatCoordStr(chr, exon_end + 1,
                                                      downstream_end)

                        (mxe_proportion1,
                         mxe_proportion2) = getMXE_MCProportion(exon_str,
                                                                upstream_start,
                                                                downstream_end,
                                                                all_jcn_count_dict,
                                                                all_coord_start2end,
                                                                all_coord_end2start)

                        # proportions used to estimate junction count for
                        # mutually exclusive events based on all junctions on
                        # each side of the exon
                        upstrm_jcn_sum2 = 0
                        dwnstrm_jcn_sum2 = 0

                        for this_start in all_coord_end2start[chr][exon_start - 1]:
                            upstrm_jcn_sum2 += all_jcn_count_dict[formatCoordStr(chr,
                                                                               this_start,
                                                                               exon_start - 1)][1]
                        for this_end in all_coord_start2end[chr][exon_end + 1]:
                            dwnstrm_jcn_sum2 += all_jcn_count_dict[formatCoordStr(chr, 
                                                                                exon_end + 1,   
                                                                                this_end)][1]
                                                                            
                        left_samp2_ct = int(round(upstrm_jcn_sum2 *
                                                  mxe_proportion2))
                        right_samp2_ct = int(round(dwnstrm_jcn_sum2 *
                                                   mxe_proportion2))

                        if norm2:
                            left_samp2_ct = int(round(left_samp2_ct/norm2))
                            right_samp2_ct = int(round(right_samp2_ct/norm2))

                        isoform_len = (2 * jcn_seq_len) + exon_end - exon_start + 1

                        left_raw_ct = left_samp2_ct
                        left_lenNorm_ct = normalizeByLen(left_samp2_ct, isoform_len)

                        right_raw_ct = right_samp2_ct
                        right_lenNorm_ct = normalizeByLen(right_samp2_ct, isoform_len)

                        raw_ct = left_raw_ct + right_raw_ct
                        lenNorm_ct = left_lenNorm_ct + right_lenNorm_ct

                        detailed_jcns.append((left_jcn_str, right_jcn_str))
                        detailed_inclusion_cts.append(((left_raw_ct, left_lenNorm_ct),
                                                       (right_raw_ct, right_lenNorm_ct)))

                    exon_cts = (raw_ct, lenNorm_ct)
                    inclusion_cts.append(exon_cts)

                    # Add to total counts
                    total_me_cts[0] += raw_ct
                    total_me_cts[1] += lenNorm_ct

                    # Add to cassette exon dict
                    if chr in cassette_exons:
                        cassette_exons[chr].add((exon_start, exon_end))
                    else:
                        cassette_exons[chr] = set([(exon_start, exon_end)])

                    # Add to excl junctions
                    if chr in excl_jcns:
                        excl_jcns[chr].add((upstream_start, exon_start - 1))
                    else:
                        excl_jcns[chr] = set([(upstream_start, exon_start - 1)])    
                        
                    excl_jcns[chr].add((exon_end + 1, downstream_end))

                    # Check for novel introns
                    if chr not in annotated_introns:
                        isNovel = True
                    elif ((upstream_start, exon_start - 1, "+") not in annotated_introns[chr]) and \
                       ((upstream_start, exon_start - 1, "-") not in annotated_introns[chr]):
                        isNovel = True
                    elif ((exon_end + 1, downstream_end, "+") not in annotated_introns[chr]) and \
                         ((exon_end + 1, downstream_end, "-") not in annotated_introns[chr]):
                        isNovel = True

                # Now print out for each exon
                for i in range(len(exon_list)):
                    # If there are only 2 exons then, only print out counts for one
                    # of them.

                    if len(exon_list) == 2 and\
                        i == 1:
                        break

                    n_or_k = "K"
                    if isNovel:
                        n_or_k = "N"

                    excl_file_raw_count = 0
                    incl_file_raw_count = 0
                    excl_file_lenNorm_count = 0
                    incl_file_lenNorm_count = 0

                    incl_file_raw_count = inclusion_cts[i][0]
                    excl_file_raw_count = total_me_cts[0] - incl_file_raw_count
                    incl_file_lenNorm_count = inclusion_cts[i][1]
                    excl_file_lenNorm_count = total_me_cts[1] - incl_file_lenNorm_count

                    # Add exon counts
    #                if full_exon_count_dict is not None:
    #                    if exon_strs[i] in full_exon_count_dict:
    #                        incl_file1_count += \
    #                                full_exon_count_dict[exon_strs[i]][0]
    #                        incl_file2_count += \
    #                                full_exon_count_dict[exon_strs[i]][1]
    #                        excl_file1_count -=\
    #                                full_exon_count_dict[exon_strs[i]][0]
    #                        excl_file2_count -=\
    #                                full_exon_count_dict[exon_strs[i]][1]

                    # Add Exclusion or Inclusion Annotation
                    # Checks percent exclusion from both files.
#                   e_or_i = checkExclusionInclusion(excl_file1_count,
#                                                   incl_file1_count,
#                                                   excl_file2_count,
#                                                   incl_file2_count) 

                    gene_name = inferGeneName(annotated_genes_by_strand, chr, strand,
                                              exon_strs)


                    out_str = "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d" %\
                                                                (n_or_k,
                                                                 "?",
                                                                 gene_name,
                                                                 chr,
                                                                 strand,
                                                                 upstream_start,
                                                                 downstream_end,
                                                                 exon_strs[i],
                                                                 ",".join(exon_strs),
                                                                 excl_file_raw_count,
                                                                 incl_file_raw_count,
                                                                 excl_file_lenNorm_count,
                                                                 incl_file_lenNorm_count)

                    if printExonCoords:
                        out_str += "\t%s\t" % exon_strs[i]

                        exon_coords.add(convertCoordStr(exon_strs[i]))                

                    out_exons = []
                    for j in range(len(exon_strs)):
                        if i == j:
                            continue
                        out_exons.append(exon_strs[j])
                        if printExonCoords:
                            exon_coords.add(convertCoordStr(exon_strs[j]))

                    if printExonCoords:
                        out_str += ",".join(out_exons)

                    me_out.write(out_str + "\n") 

                    # Write to all event file
                    excl_jcn_str_list = []
                    excl_jcn_ct_strs_raw = []
                    for j in range(len(detailed_jcns)):
                        if i == j:
                            continue
                        jcns_tuple = detailed_jcns[j]
                        excl_jcn_str_list.append(",".join(jcns_tuple))

                        left_jcn_ct_raw = detailed_inclusion_cts[j][0][0]
                        right_jcn_ct_raw = detailed_inclusion_cts[j][1][0]

                        excl_jcn_ct_strs_raw.append("%d,%d" % (left_jcn_ct_raw, right_jcn_ct_raw))
                        

                    out_str = getAllEventStr(n_or_k,
                                             "mutually_exclusive", 
                                             gene_name, chr, strand, 
                                             ";".join(excl_jcn_str_list),
                                             ",".join(detailed_jcns[i]),
                                             ";".join(out_exons),
                                             exon_strs[i],
                                             "",
                                             ";".join(const_strs),
                                             ";".join(excl_jcn_ct_strs_raw),
                                             "%d,%d" % (detailed_inclusion_cts[i][0][0],detailed_inclusion_cts[i][1][0]),
                                             repr(excl_file_raw_count),
                                             repr(incl_file_raw_count),
                                             "","",
                                             None, None,
                                             "",
                                             None,
                                             "",
                                             None)
                                             

                    all_event_info_out.write(out_str + "\n")


                                                             
def removeOverlappingCassetteExons(cassette_exons, start_or_end):
    """
    Iterates through every exon in the set and checks for an groups of
    overlapping exons.  Finally, all overlapping exons in the set are removed.

    If the cassette exons have the same start_or_end position, then they are
    not removed.
    """
    for chr in cassette_exons:
        exons2remove = set([])
        exon_list = list(cassette_exons[chr])
        exon_list.sort()

        for i in range(len(exon_list)-1):
            if exon_list[i] in exons2remove:
                continue
            this_overlap_set = set([exon_list[i]])

            for j in range(i+1, len(exon_list)): 
                # check if this coord overlaps with any existing coord in
                # this_overlap set
                if coordOverlaps_wSet(exon_list[j], this_overlap_set):
                    this_overlap_set.add(exon_list[j])

                else: 
                    # Since list is sorted, the overlapping exons should be
                    # adjacent to each other.
                    break

            if len(this_overlap_set) > 1:
                if not allHaveSameCoord(this_overlap_set, start_or_end):
                    # Full overlap set was made, now add to exons2remove set
                    exons2remove = exons2remove.union(this_overlap_set)

        # Remove the exons
        for coord in exons2remove:
            cassette_exons[chr].remove(coord)

def removeRedundantSets(coord_set):
    """ 
    Removes any coord_sets that are subset of other coord_sets
    """
    # Remove singletons first
    new_set = []
    for elem in coord_set:
        if len(elem) != 1:
            new_set.append(elem)
       
    non_redundant_set = list(new_set)
    for i in range(len(new_set)):
        for j in range(len(new_set)):
            if i == j:
                continue
            if new_set[i].issubset(new_set[j]):
                non_redundant_set[i] = None

    return_set = []
    for elem in non_redundant_set:
        if elem is not None:
            return_set.append(elem)

    return return_set

def splitSharedRegions(alt_start_or_end, 
                       cts1_list, cts2_list,
                       proportions1, proportions2):
    """
    Returns updated counts but to each isoform: isoform_cts1,
                                                isoform_cts2
    
    The shared region is determined by whether it is an alternative start or
    end splice site and the proportions given
    """
    if len(cts1_list) != len(cts2_list):
        ERROR_LOG.write("Something wrong in splitSharedRegions.\n")
        return ([0]*len(cts1_list), [0]*len(cts2_list))
    if len(proportions1) != len(proportions2):
        ERROR_LOG.write("Something wrong in splitSharedRegions.\n")
        return ([0]*len(cts1_list), [0]*len(cts2_list))
    if len(cts1_list) != len(proportions2):
        ERROR_LOG.write("Something wrong in splitSharedRegions.\n")
        return ([0]*len(cts1_list), [0]*len(cts2_list))

    num_isoforms = len(cts1_list)

    if num_isoforms == 1:
        return list(cts1_list), list(cts2_list)

    isoform_cts1 = [0] * num_isoforms
    isoform_cts2 = [0] * num_isoforms

    if alt_start_or_end == "alt_start":
        # First count is shared by all isoforms.
        for cts_pos in range(num_isoforms):
            propor_subset1 = proportions1[cts_pos:]
            propor_subset2 = proportions2[cts_pos:]

            updated_proportions1 = [0] * cts_pos
            for p1 in proportions1[cts_pos:]:
                if sum(propor_subset1) == 0:
                    updated_proportions1.append(0)
                else:
                    updated_proportions1.append(p1/sum(propor_subset1))
            updated_proportions2 = [0] * cts_pos
            for p2 in proportions2[cts_pos:]:
                if sum(propor_subset2) == 0:
                    updated_proportions2.append(0)
                else:
                    updated_proportions2.append(p2/sum(propor_subset2))
            
            for isoform_pos in range(cts_pos, num_isoforms):
                # Add relative proportion of counts to the isoform counts
                isoform_cts1[isoform_pos] += int(round(cts1_list[cts_pos] *
                                                   updated_proportions1[isoform_pos]))
                isoform_cts2[isoform_pos] += int(round(cts2_list[cts_pos] *
                                                   updated_proportions2[isoform_pos]))
                
        
    else: # alt_end
        # Last count is shared by all isoforms
        for cts_pos in range(num_isoforms): 
            propor_subset1 = proportions1[:cts_pos+1]
            propor_subset2 = proportions2[:cts_pos+1]

            updated_proportions1 = []
            for p1 in proportions1[:cts_pos+1]:
                if sum(propor_subset1) == 0:
                    updated_proportions1.append(0)
                else:
                    updated_proportions1.append(p1/sum(propor_subset1))
            updated_proportions2 = []
            for p2 in proportions2[:cts_pos+1]:
                if sum(propor_subset2) == 0:
                    updated_proportions2.append(0)
                else:
                    updated_proportions2.append(p2/sum(propor_subset2))
    
            for isoform_pos in range(cts_pos + 1):
                # Add relative proportion of counts to isoform counts
                isoform_cts1[isoform_pos] += int(round(cts1_list[cts_pos] *
                                                   updated_proportions1[isoform_pos]))
                isoform_cts2[isoform_pos] += int(round(cts2_list[cts_pos] *
                                                   updated_proportions2[isoform_pos]))

    return isoform_cts1, isoform_cts2


def subtract_vectors(vec1, vec2):
    if len(vec1) != len(vec2):
        ERROR_LOG.write("Problem in subtract_vectors. Vectors are not the same length.\n")
        return None

    new_vec = []
    for i in range(len(vec1)):
        new_vec.append(vec1[i] - vec2[i])

    return new_vec

def sumExclusion_Inclusion_counts(file_str, 
                                  ce2total_counts,
                                  alt_donor2total_counts,
                                  alt_accept2total_counts,
                                  afe2total_counts,
                                  ale2total_counts,
                                  mxe2total_counts,
                                  mc2total_counts,
                                  printExonCoords,
                                  norm1, norm2, jcn_seq_len):
    """
    Also serves as a check point.

    For cassette exon events, counts come from previous calculations

    For alternative donor an acceptor events, the counts come from previous
    calculations.
    """
    file = open(file_str)
    lines = file.readlines()
    file.close() 

    file2 = open(file_str, "w")

    for line in lines:
        line = line.rstrip("\n")

        line_elems = line.split("\t") 

        type = line_elems[1]
    
        # Check for required counts
        required_indices = []
        if type == "cassette":
            required_indices = [13,14]
        elif type == "mutually_exclusive":
            required_indices = [13,14]
        elif type == "coord_cassette":
            required_indices = [13,14]
        elif type == "alternative_donor":
            required_indices = [13,14]
        elif type == "alternative_acceptor":
            required_indices = [13,14]
        elif type == "jcn_only_AD":
            required_indices = [13,14]
        elif type == "jcn_only_AA":
            required_indices = [13,14]
        elif type == "alternative_first_exon":
            required_indices = [13,14]
        elif type == "alternative_last_exon":
            required_indices = [13,14]
        elif type == "intron_retention":
            required_indices = [13, 20]
        else:
            print "Unknown type in: %s, %s" % (type, file_str)

        if printExonCoords:
            if type == "cassette":
                required_indices += [18]
            elif type == "mutually_exclusive":
                required_indices += [17,18]
            elif type == "coord_cassette":
                required_indices += [18]
#           elif type == "alternative_donor":
#               required_indices += [25, 27]
#           elif type == "alternative_acceptor":
#               required_indices += [25, 27]
            
        
        for required_idx in required_indices:
            if line_elems[required_idx] == "":
                print "Problem with %s" % file_str
                print "Line: %s" % line
                sys.exit(1)

        sum_excl_raw = 0
        sum_incl_raw = 0
        sum_excl_lenNorm = 0
        sum_incl_lenNorm = 0
            
   
        if printExonCoords:
            excl_raw_cols = [13, 17] 

            incl_raw_cols = [14, 18, 20]

        else:
            excl_raw_cols = [13] 

            incl_raw_cols = [14, 20]

        if type == "cassette":
            event_key = line_elems[8]

            if event_key not in ce2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find CE key. %s\n" % event_key)

            sum_excl_raw = ce2total_counts[event_key][0]
            sum_incl_raw = ce2total_counts[event_key][1]
            sum_excl_lenNorm = ce2total_counts[event_key][2]
            sum_incl_lenNorm = ce2total_counts[event_key][3]

        elif type == "alternative_donor":
            event_key = (line_elems[6], line_elems[5])
            if event_key not in alt_donor2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find AD key. %s\n" % event_key)

            sum_excl_raw = alt_donor2total_counts[event_key][0]
            sum_incl_raw = alt_donor2total_counts[event_key][1]
            sum_excl_lenNorm = alt_donor2total_counts[event_key][2]
            sum_incl_lenNorm = alt_donor2total_counts[event_key][3]
        elif type == "alternative_acceptor":
            event_key = (line_elems[6], line_elems[5])
            if event_key not in alt_accept2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find AA key. %s\n" % event_key)

            sum_excl_raw = alt_accept2total_counts[event_key][0]
            sum_incl_raw = alt_accept2total_counts[event_key][1]
            sum_excl_lenNorm = alt_accept2total_counts[event_key][2]
            sum_incl_lenNorm = alt_accept2total_counts[event_key][3]
        elif type == "jcn_only_AD" or type == "jcn_only_AA":
            # Only junction counts are used for the final quantification
            sum_excl_raw = int(line_elems[13])
            sum_incl_raw = int(line_elems[14])
            sum_excl_lenNorm = int(round((float(sum_excl_raw))/(jcn_seq_len/DEF_EXON_LEN_NORM)))
            sum_incl_lenNorm = int(round((float(sum_incl_raw))/(jcn_seq_len/DEF_EXON_LEN_NORM)))
        elif type == "alternative_first_exon":
            event_key = (line_elems[6], line_elems[5])
            if event_key not in afe2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find AFE key. %s\n" % event_key)

            sum_excl_raw = afe2total_counts[event_key][0]
            sum_incl_raw = afe2total_counts[event_key][1]
            sum_excl_lenNorm = afe2total_counts[event_key][2]
            sum_incl_lenNorm = afe2total_counts[event_key][3]

        elif type == "alternative_last_exon":
            event_key = (line_elems[6], line_elems[5])
            if event_key not in ale2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find ALE key. %s\n" % event_key)

            sum_excl_raw = ale2total_counts[event_key][0]
            sum_incl_raw = ale2total_counts[event_key][1]
            sum_excl_lenNorm = ale2total_counts[event_key][2]
            sum_incl_lenNorm = ale2total_counts[event_key][3]

        elif type == "mutually_exclusive":
            event_key = (line_elems[8], line_elems[7])
            if event_key not in mxe2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find MXE key. %s\n" % event_key)

            sum_excl_raw = mxe2total_counts[event_key][0]
            sum_incl_raw = mxe2total_counts[event_key][1]
            sum_excl_lenNorm = mxe2total_counts[event_key][2]
            sum_incl_lenNorm = mxe2total_counts[event_key][3]
        elif type == "coord_cassette":
            event_key = (line_elems[8], line_elems[5])
            if event_key not in mc2total_counts:
                ERROR_LOG.write("sumExclusion_Inclusion_counts: cannot find MC key. %s\n" % event_key)

            sum_excl_raw = mc2total_counts[event_key][0]
            sum_incl_raw = mc2total_counts[event_key][1]
            sum_excl_lenNorm = mc2total_counts[event_key][2]
            sum_incl_lenNorm = mc2total_counts[event_key][3]
        else: # intron retention
            sum_excl_raw = getColSums(excl_raw_cols, line_elems)
            sum_incl_raw = getColSums(incl_raw_cols, line_elems)

            sum_excl_lenNorm = normalizeByLen(sum_excl_raw, jcn_seq_len)
            sum_incl_lenNorm = normalizeByLen(sum_incl_raw, 2 * jcn_seq_len)

        if hasNegativeVals(sum_excl_raw, 
                           sum_incl_raw, 
                           sum_excl_lenNorm,
                           sum_incl_lenNorm):
            # ERROR WAS ALREADY PRINTED
            sum_excl_raw = 0
            sum_incl_raw = 0
            sum_excl_lenNorm = 0
            sum_incl_lenNorm = 0

        line_elems.append(repr(sum_excl_raw))
        line_elems.append(repr(sum_incl_raw))
        line_elems.append(repr(sum_excl_lenNorm))
        line_elems.append(repr(sum_incl_lenNorm))

        outline = "\t".join(line_elems)

        file2.write(outline + "\n")        

    file2.close()

def translateInput(bed_line):
    """
    Junction counts should be in BED format.
    Output is a coord_str chr_intron.start_intron.end
    strand
    count

    Coordinates are converted to 1-based.
    """
    input_list = bed_line.split("\t")

    # Input checking
    if len(input_list) != 12:
        print "Error with BED file format."
        sys.exit(1)

    if int(input_list[-3]) != 2:
        print "Cannot process reads that span more than one junciton."
        return None, None, None

    chr = input_list[0]

    if not chr.startswith("chr"):
        chr = "chr" + chr
#   if "" in chr:
#       chr = chr.replace("","-")

    chromStart = int(input_list[1])

    blockSizes = map(int, input_list[-2].split(","))
    blockStarts = map(int, input_list[-1].split(","))
    
    intron_start = chromStart + blockSizes[0] + 1
    intron_end = chromStart + blockStarts[1]

    jcn_str = formatCoordStr(chr, intron_start, intron_end)

    strand = input_list[5]

    if strand != "+" and strand != "-" and strand != ".":
        strand = None
        print "Something wrong with strand of junction: %s\n" % bed_line

#   if strand == ".":
#       strand = None

    count = int(input_list[4])

#    if input_list[1] == "+":
#        strand = "+"
#    if input_list[1] == "-":
#        strand = "-"

    return jcn_str, strand, count

def updateProportions(alt_start_or_end, proportions1, proportions2):
    """
    Returns updated proportions where the first or last element is not
    considered depending on whether it is an alt_start or alt_end.
    Retunrs an updated list of proportions. updated_proportions1,
    updated_proportions2
    """
    updated_proportions1 = []
    updated_proportions2 = []

    if alt_start_or_end == "alt_start":
        propor_subset1 = proportions1[1:]
        propor_subset2 = proportions2[1:]
    else:
        propor_subset1 = proportions1[:-1]
        propor_subset2 = proportions2[:-1]

    for p1 in propor_subset1:
        if sum(propor_subset1) == 0:
            updated_proportions1.append(0)
        else:
            updated_proportions1.append(p1/sum(propor_subset1))
    for p2 in propor_subset2:
        if sum(propor_subset2) == 0:
            updated_proportions2.append(0)
        else:
            updated_proportions2.append(p2/sum(propor_subset2))

    return updated_proportions1, updated_proportions2

def updateStrand(cur_strand, new_strand):
    if cur_strand == "." or new_strand == ".":
        return "."

    if not new_strand:
        return "."

    if cur_strand:
        if cur_strand == new_strand:
            return cur_strand
        else:
            return "."
    else:
        return new_strand

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
#if __name__ == "__main__": profile.run('main()')
