#!/lab/64/bin/python
# run_getASEventReadCounts_multiSample.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.
"""Runs each sample against a pseudo reference.
   Difference from previous version is allowing for different input options
   where the two different databases are given.
"""

import sys
import optparse 
import os

from subprocess import Popen
from helperFunctions import runCmd
from broad_helperFunctions import runLSF

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
                          help="Comma separated list of samples to run.",
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
    opt_parser.add_option("--method",
                          dest="method",
                          type="string",
                          help="""Type of correction method:
                                  'BH' - Benjamini & Hochberg,
                                  'bonferroni'""",
                          default=None)
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
    opt_parser.add_option("--lengthNorm",
                          dest="lengthNorm",
                          action="store_true",
                          help="""Default is to normalize read counts by
                                 isoform length. This will option will specify
                                 to not normalize by isoform length.""",
                          default=False)
    opt_parser.add_option("--fasta",
                           dest="genome_file",
                           type="string",
                           help="""Contains the genome sequence organized by
                                   chromosome.""",
                           default=None)
    opt_parser.add_option("-p",
                          dest="num_processes",
                          type="int",
                          help="""Will run getASEventReadCounts.py
                                  simultaneously with this many samples.
                                  Default=%d""" % DEF_NUM_PROCESSES,
                          default=DEF_NUM_PROCESSES)
    opt_parser.add_option("--LSF",
                          dest="run_lsf",
                          action="store_true",
                          help="""Will launch jobs on LSF. Default is running on
                                  local.""",
                          default=False)
    opt_parser.add_option("--by_chr",
                          dest="by_chr",
                          action="store_true",
                          help="""Indicates that input files are broken up by
                                  chromosome""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-s")
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")
    opt_parser.check_required("--txt_db1")
    opt_parser.check_required("--txt_db2")
    opt_parser.check_required("--method")
    opt_parser.check_required("--fasta")
    opt_parser.check_required("--jcn_seq_len")

    samples = options.samples.split(",")

    if os.path.exists(options.input_dir):
        input_dir = os.path.abspath(options.input_dir)
    else:
        print "Input directory does not exist."
        opt_parser.print_help()
        sys.exit(1)

    if input_dir.endswith("/"):
        input_dir = input_dir.rstrip("/")

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

    method = options.method
    genome_file = os.path.abspath(options.genome_file)

    jcn_seq_len = options.jcn_seq_len

    num_processes = options.num_processes
    run_LSF = options.run_lsf

    by_chr = options.by_chr

    if by_chr:
        chr_list = getChr(input_dir)       
        
        ctr = 0
        for samp in samples:
            ctr += 1
            # Check for output subdirectory
            samp_dir = output_dir + "/" + samp
            if not os.path.exists(samp_dir):
                os.mkdir(samp_dir)

            for chr in chr_list:
                chr_dir = samp_dir + "/" + samp + "_" + chr
                if not os.path.exists(chr_dir):
                    os.mkdir(chr_dir)

                os.chdir(chr_dir)

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
                cmd += "--method %s " % method
                cmd += "--jcn_seq_len %d " % jcn_seq_len
                cmd += "--fasta %s " % genome_file
                cmd += "--by_chr %s " % chr

                if options.lengthNorm:
                    cmd += "--lengthNorm "

                # Now for databases
                if options.sqlite_db_dir:
                    cmd += "--sqlite_db_dir %s" % options.sqlite_db_dir
                else: # use MySQL
                    if options.passwd == "":
                        cmd += "--host %s --user %s" % (options.host,
                                                        options.user)
                    else:
                        cmd += "--host %s --user %s --passwd %s" % (options.host,
                                                                 options.user,
                                                                 options.passwd)

                if run_LSF:
                    runLSF(cmd, 
                           "%s_%s.getASEventReadCounts.bsub.out" % (samp, chr),
                           samp,
                           "hour") 
                    continue

                if ctr % num_processes == 0:
                    os.system(cmd)
                else:
                    print cmd
                    Popen(cmd, shell=True, executable=SHELL)

    else:
        ctr = 0
        for samp in samples:
            ctr += 1
            # Check for output subdirectory
            full_output_dir = output_dir + "/" + samp
            if not os.path.exists(full_output_dir):
                os.mkdir(full_output_dir)

            os.chdir(full_output_dir)

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
            cmd += "--method %s " % method
            cmd += "--jcn_seq_len %d " % jcn_seq_len
            cmd += "--fasta %s " % genome_file

            if options.lengthNorm:
                cmd += "--lengthNorm "

            # Now for databases
            if options.sqlite_db_dir:
                cmd += "--sqlite_db_dir %s" % options.sqlite_db_dir
            else: # use MySQL
                if options.passwd == "":
                    cmd += "--host %s --user %s" % (options.host,
                                                    options.user)
                else:
                    cmd += "--host %s --user %s --passwd %s" % (options.host,
                                                             options.user,
                                                             options.passwd)

            if run_LSF:
                runLSF(cmd, 
                       "%s.getASEventReadCounts.bsub.out" % samp,
                       samp,
                       "week") # Week cue if running whole samples
                continue

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
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
