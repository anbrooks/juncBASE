#!/lab/64/bin/python
# run_getASEventReadCounts_multiSample.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.
"""Runs each sample against a pseudo reference.
   Difference from previous version is allowing for different input options
m
   where the two different databases are given.
"""

import sys
import optparse 
import os
import pdb

from subprocess import Popen
from helperFunctions import runCmd, runLSF

from createPseudoSample import getChr
#############
# CONSTANTS #
#############
BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/getASEventReadCounts.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: getASEventReadCounts.py needs to be in the same directory."
    sys.exit(1)

SHELL = "/bin/tcsh"

DEF_NUM_PROCESSES = 2

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
    opt_parser.add_option("-s",
                          dest="samples",
                          type="string",
                          help="""Comma separated list of samples to run or a file
                                  with the names of the samples to run. If a
                                  file is given, the first column will be used
                                  as the sample column and is assumed
                                  tab-delimited""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="input_dir",
                          type="string",
                          help="""Root of directory containing input files for all
                                  samples.""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_dir",
                          type="string",
                          help="""Root of directory that will place all output
                                  into subdirectories.""",
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
                                  names and whether an intron/junction is
                                  annotated or not. By default, txt_db1 will be used for this
                                  information.""",
                          default=None)
#   opt_parser.add_option("--method",
#                         dest="method",
#                         type="string",
#                         help="""Type of correction method:
#                                 'BH' - Benjamini & Hochberg,
#                                 'bonferroni'""",
#                         default=None)
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
#   opt_parser.add_option("--lengthNorm",
#                         dest="lengthNorm",
#                         action="store_true",
#                         help="""Default is to normalize read counts by
#                                isoform length. This will option will specify
#                                to not normalize by isoform length.""",
#                         default=False)
#   opt_parser.add_option("--fasta",
#                          dest="genome_file",
#                          type="string",
#                          help="""Contains the genome sequence organized by
#                                  chromosome.""",
#                          default=None)
    opt_parser.add_option("-p",
                          dest="num_processes",
                          type="int",
                          help="""Will run getASEventReadCounts.py
                                  simultaneously with this many samples.
                                  Default=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)
    opt_parser.add_option("--nice",
                          dest="nice",
                          action="store_true",
                          help="When running locally, use nice",
                         default=False)
    opt_parser.add_option("--LSF",
                          dest="run_lsf",
                          action="store_true",
                          help="""Will launch jobs on LSF. Default is running on
                                  local.""",
                          default=False)
    opt_parser.add_option("--week",
                          dest="week",
                          action="store_true",
                          help="Will launch jobs on LSF using week queue.",
                          default=False)
    opt_parser.add_option("--by_chr",
                          dest="by_chr",
                          action="store_true",
                          help="""Indicates that input files are broken up by
                                  chromosome""",
                          default=False)
    opt_parser.add_option("--force",
                          dest="force",
                          action="store_true",
                          help="""By default, will check for the existence of
                                  the final output before running commands. This
                                  option will force all runs.""",
                          default=False)
    opt_parser.add_option("--check",
                          dest="check",
                          action="store_true",
                          help="""Will check samples that are not done and print
                                  out which need to still be run""",
                         default=False)
    opt_parser.add_option("--print_cmd",
                          dest="print_cmd",
                          action="store_true",
                          help="""Will print commands that will be run, but will
                                  not run them. Used for debugging.""",
                         default=False)
    opt_parser.add_option("--keep_intermediate",
                           dest="keep_interm",
                           action="store_true",
                           help="""Will remove intermediate files by default.
                                   Use this option to keep them.""",
                           default=False)


    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-s")
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")
    opt_parser.check_required("--txt_db1")
    opt_parser.check_required("--txt_db2")
#    opt_parser.check_required("--method")
#    opt_parser.check_required("--fasta")
    opt_parser.check_required("--jcn_seq_len")

    samples = getSampleNames(options.samples)

    if os.path.exists(options.input_dir):
        input_dir = os.path.abspath(options.input_dir)
    else:
        print "Input directory does not exist."
        opt_parser.print_help()
        sys.exit(1)

    if input_dir.endswith("/"):
        input_dir = input_dir.rstrip("/")

    if options.sqlite_db_dir:
        sqlite_db_dir = formatDir(options.sqlite_db_dir)

    if os.path.exists(options.output_dir):
        output_dir = os.path.abspath(options.output_dir)
    else:
        os.mkdir(options.output_dir)
        output_dir = os.path.abspath(options.output_dir)
        print "Creating output directory: %s" % output_dir

    if output_dir.endswith("/"):
        output_dir = output_dir.rstrip("/")


    txt_db1 = options.txt_db1
    txt_db2 = options.txt_db2
    txt_db3 = options.txt_db3

#    method = options.method
#    genome_file = os.path.abspath(options.genome_file)

    jcn_seq_len = options.jcn_seq_len

    num_processes = options.num_processes
    run_LSF = options.run_lsf

    keep_interm = options.keep_interm

    week = options.week

    nice = options.nice

    print_cmd = options.print_cmd 

    force = options.force
    check = options.check

    by_chr = options.by_chr

    if by_chr:
        chr_list = getChr(input_dir)       
        
        ctr = 0
        for samp in samples:
            # Check for output subdirectory
            samp_dir = output_dir + "/" + samp
            if not os.path.exists(samp_dir):
                os.mkdir(samp_dir)

            for chr in chr_list:
                chr_dir = samp_dir + "/" + samp + "_" + chr
                if not os.path.exists(chr_dir):
                    os.mkdir(chr_dir)

                os.chdir(chr_dir)

                expected_out_file = "%s_%s_finished.txt" % (samp, chr)

                file_is_present = False

                try:
                    if os.path.getsize(expected_out_file) != 0:
                        file_is_present = True
                    else:
                        if check:
                            print "File is empty: %s,%s,%s" % (samp, chr, expected_out_file)
                except:                                                              
                    if check:
                        print "Does not exist: %s,%s,%s" % (samp, chr,
                                                            expected_out_file)

                if check:
                    continue

                if force:
                    # Delete previous expected out file to prevent confusions
                    # with previous runs when rechecking again
                    if file_is_present:
                        os.system("rm " + expected_out_file)
                else:
                    if file_is_present:
                        continue

                ctr += 1

                cmd = "python %s " % SCRIPT
                cmd += "--jcn1 %s/pseudo_%s/pseudo_%s_junctions.bed " % (input_dir,
                                                                         chr,
                                                                         chr)
                cmd += "--jcn2 %s/%s/%s_%s/%s_%s_junctions.bed " % (input_dir, 
                                                                    samp, 
                                                                    samp, chr,
                                                                    samp, chr)
                cmd += "--genome_reads1 %s/pseudo_%s/pseudo_%s_genome_reads.txt.gz " % (input_dir,
                                                                                        chr,
                                                                                        chr)
                cmd += "--genome_reads2 %s/%s/%s_%s/%s_%s_genome_reads.txt.gz " % (input_dir, 
                                                                                   samp, 
                                                                                   samp, chr,
                                                                                   samp, chr)
                cmd += "--ie1 %s/pseudo_%s/pseudo_%s_intron_exon_junction_counts.txt " % (input_dir,
                                                                                          chr,
                                                                                          chr) 
                cmd += "--ie2 %s/%s/%s_%s/%s_%s_intron_exon_junction_counts.txt " % (input_dir, 
                                                                                     samp, 
                                                                                     samp, chr,
                                                                                     samp, chr)
                cmd += "-p %s_%s " % (samp, chr)
                cmd += "--txt_db1 %s " % txt_db1
                cmd += "--txt_db2 %s " % txt_db2
                if txt_db3:
                    cmd += "--txt_db3 %s " % txt_db3
            
#                cmd += "--method %s " % method
                cmd += "--jcn_seq_len %d " % jcn_seq_len
#                cmd += "--fasta %s " % genome_file
                cmd += "--by_chr %s " % chr

                if keep_interm:
                    cmd += "--keep_intermediate "

                # Now for databases
                if options.sqlite_db_dir:
                    cmd += "--sqlite_db_dir %s" % sqlite_db_dir
                else: # use MySQL
                    if options.passwd == "":
                        cmd += "--host %s --user %s" % (options.host,
                                                        options.user)
                    else:
                        cmd += "--host %s --user %s --passwd %s" % (options.host,
                                                                 options.user,
                                                                 options.passwd)
                if print_cmd:
                    if not run_LSF:
                        if nice:
                            cmd = "nice " + cmd
                    print cmd
                    continue

                if run_LSF:
                    if week:
                        queue = "week"
                    else:
                        queue = "hour"

                    runLSF(cmd, 
                           "%s_%s.getASEventReadCounts.bsub.out" % (samp, chr),
                           samp + "_" + chr,
                           queue) 
                    continue

                if nice:
                    cmd = "nice " + cmd

                if ctr % num_processes == 0:
                    os.system(cmd)
                else:
                    print cmd
                    Popen(cmd, shell=True, executable=SHELL)

    else:
        ctr = 0
        for samp in samples:
            # Check for output subdirectory
            full_output_dir = output_dir + "/" + samp
            if not os.path.exists(full_output_dir):
                os.mkdir(full_output_dir)

            os.chdir(full_output_dir)

            expected_out_file = "%s_finished.txt" % samp
                                                                               
            file_is_present = False                                            

            try:
                if os.path.getsize(expected_out_file) != 0:
                    file_is_present = True
                else:
                    if check:
                        print "File is empty: %s,%s" % (samp, expected_out_file)
            except:                                                              
                if check:
                    print "Does not exist: %s, %s" % (samp, expected_out_file)
            
            if check:
                continue
            
            if force:
                if file_is_present:
                    os.system("rm " + expected_out_file)
            else:
                if file_is_present:
                    continue

            ctr += 1

            cmd = "python %s " % SCRIPT
            cmd += "--jcn1 %s/pseudo/pseudo_junctions.bed " % input_dir
            cmd += "--jcn2 %s/%s/%s_junctions.bed " % (input_dir, samp, samp)
            cmd += "--genome_reads1 %s/pseudo/pseudo_genome_reads.txt.gz " % input_dir
            cmd += "--genome_reads2 %s/%s/%s_genome_reads.txt.gz " % (input_dir, samp, samp)
            cmd += "--ie1 %s/pseudo/pseudo_intron_exon_junction_counts.txt " % input_dir
            cmd += "--ie2 %s/%s/%s_intron_exon_junction_counts.txt " % (input_dir, samp, samp)
            cmd += "-p %s " % samp
            cmd += "--txt_db1 %s " % txt_db1
            cmd += "--txt_db2 %s " % txt_db2
            if txt_db3:
                cmd += "--txt_db3 %s " % txt_db3

#            cmd += "--method %s " % method
            cmd += "--jcn_seq_len %d " % jcn_seq_len
#            cmd += "--fasta %s " % genome_file

            if keep_interm:
                cmd += "--keep_intermediate "

            # Now for databases
            if options.sqlite_db_dir:
                cmd += "--sqlite_db_dir %s" % sqlite_db_dir
            else: # use MySQL
                if options.passwd == "":
                    cmd += "--host %s --user %s" % (options.host,
                                                    options.user)
                else:
                    cmd += "--host %s --user %s --passwd %s" % (options.host,
                                                             options.user,
                                                             options.passwd)

            if print_cmd:
                if not run_LSF:
                    if nice:
                        cmd = "nice " + cmd
                print cmd
                continue

            if run_LSF:
                runLSF(cmd, 
                       "%s.getASEventReadCounts.bsub.out" % samp,
                       samp,
                       "week") # Week cue if running whole samples
                continue

            if nice:
                cmd = "nice " + cmd

            if ctr % num_processes == 0:
                os.system(cmd)
            else:
                print cmd
                Popen(cmd, shell=True, executable=SHELL)
            
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getSampleNames(samples_option):
    if not os.path.exists(samples_option):
#        if "," in samples_option:
        return samples_option.split(",")
#        print "Cannot find sample names file: %s" % samples_option
#        sys.exit(1)

    s_file = open(samples_option)


    samples = []
    for line in s_file:
        line = formatLine(line)
        lineList = line.split("\t")
        samples.append(lineList[0]) 

    return samples
    
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
