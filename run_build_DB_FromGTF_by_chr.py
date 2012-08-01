#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# run_build_DB_FromGTF_by_chr.py
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

from subprocess import Popen

from pysqlite_wrap import DB
from helperFunctions import updateDictOfLists, runCmd
from broad_helperFunctions import runLSF

#############
# CONSTANTS #
#############
# Database variables
DB_DIR = "./"

BIN_DIR = os.path.realpath(os.path.dirname(sys.argv[0]))
SCRIPT = "%s/build_DB_FromGTF.py" % BIN_DIR
if not os.path.exists(SCRIPT):
    print "ERROR: build_DB_FromGTF.py needs to be in the same directory."
    sys.exit(1)

SHELL = "/bin/tcsh"

DEF_NUM_PROCESSES = 2
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
    opt_parser.add_option("--initialize",
                          dest="initialize",
                          action="store_true",
                          help="""Will split up the gtf file into separate temp files
                                  and initalize the database.""",
                          default=False)
    opt_parser.add_option("--tmp_dir",
                          dest="tmp_dir",
                          type="string",
                          help="""Directory to place temporary files and to look
                                  for temporary files.""",
                          default=None)
    opt_parser.add_option("--keep_temp",
                          dest="keep_temp",
                          action="store_true",
                          help="""TEMP FILES ARE KEPT FOR NOW. Will keep the temporary gtf files. Default is
                                  to delete them.""",
                          default=False)
    opt_parser.add_option("-g",
                          dest="gtf_file",
                          type="string",
                          help="GTF annotation file.",
                          default=None)
    opt_parser.add_option("--use_gene_name",
                          dest="use_gene_name",
                          action="store_true",
                          help="""By default, the gene_id attribute will be used
                                  for the gene name used in the database, but
                                  the gene_name attribute can be used
                                  instead.""",
                          default=False)
    # May revisit this option, but do not need now
#   opt_parser.add_option("-f",
#                         dest="genome_file_name",
#                         type="string",
#                         help="""Fasta file containing all chromosome
#                                 sequences.  If this option is given, exon and
#                                 intron sequences will be stored in the
#                                 database as well. Chromosome names must be the
#                                 same format as in the gtf file.""",
#                         default=None)
    opt_parser.add_option("-d",
                          dest="db_name",
                          type="string",
                          help="Name of the new database",
                          default=None)
    opt_parser.add_option("--sqlite_db_dir",
                          dest="sqlite_db_dir",
                          type="string",
                          help="Location to put sqlite database. Default=%s" % DB_DIR,
                          default=DB_DIR)
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


    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-g")
    opt_parser.check_required("--tmp_dir")
    opt_parser.check_required("-d")

    gtf_file_name = options.gtf_file
    tmp_dir = formatDir(options.tmp_dir)

    db_name = options.db_name

    sqlite_db_dir = options.sqlite_db_dir

    num_processes = options.num_processes
    run_lsf = options.run_lsf

    force = options.force
    check = options.check
    print_cmd = options.print_cmd

    ##############
    # INITIALIZE #
    ##############

    # If it's initilalizing, split gtf file and initialize database return
    if options.initialize:    
        chr2lines = {}

        gtf_file_path = gtf_file_name
        gtf_file_name = gtf_file_name.split("/")[-1]
        gtf_file_comp = gtf_file_name.split(".")   
        gtf_file_prefix = ".".join(gtf_file_comp[:-1])
 
        gtf_file = open(gtf_file_path)

        for line in gtf_file:
            this_chr = line.split("\t")[0]
            updateDictOfLists(chr2lines, this_chr, line)
        gtf_file.close()

        for chr in chr2lines:
            tmp_chr_file = open("%s/%s_%s.gtf" % (tmp_dir,
                                                  gtf_file_prefix, chr),
                                "w")
            for line in chr2lines[chr]:
                tmp_chr_file.write(line)
            tmp_chr_file.close()

        # Now initialize the database
        cmd = "python %s " % SCRIPT
        cmd += "--initialize -d %s" % db_name
        os.system(cmd)
        
        sys.exit(0)

    ##################
    # BUILD DATABASE #
    ##################
    db = DB(sqlite_db_dir)

    # Use gtf file to figure out temp file names, Build the database from them
    tmp_file_list = []
    
    gtf_file_name = gtf_file_name.split("/")[-1]
    gtf_file_comp = gtf_file_name.split(".")   
    gtf_file_prefix = ".".join(gtf_file_comp[:-1])

    for this_file in os.listdir(tmp_dir):
        if gtf_file_prefix in this_file:
            if this_file == gtf_file_name:
                continue
            tmp_file_list.append(this_file)


    # Now run script for every chromosome file
    ctr = 0
    for tmp_file in tmp_file_list:

        this_chr = getChr(tmp_dir + "/" + tmp_file)

        if (not force) or check:

            # For now, just checks that records exist in the database, It is
            # better to force since it difficult to really know if a chromosome was
            # built or not.
            chr_built = checkChr(db, db_name, this_chr)
            
            if chr_built:
                if not force:
                    continue

            if check:
                if not chr_built:
                    print "Chromosome %s not built" % this_chr
                    continue

        ctr += 1

        cmd = "python %s " % SCRIPT
        cmd += "-g %s/%s " % (tmp_dir, tmp_file)
        cmd += "-d %s " % db_name

        if options.use_gene_name:
            cmd += "--use_gene_name "

        cmd += "--sqlite_db_dir %s" % sqlite_db_dir
        
        if print_cmd:
            print cmd
            continue


        if run_lsf:
            runLSF(cmd,
                   "%s.build_DB.bsub.out" % this_chr,
                   this_chr + "build_DB",
                   "hour")
            continue
        
        if ctr % num_processes == 0:                             
            os.system(cmd)                                       
        else:
            print cmd                                       
            Popen(cmd, shell=True, executable=SHELL)     
    
    # Remove temp files, but first check that exons are returned from the same
    # chromosome in the database
#    if not options.keep_temp:
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def checkChr(db, db_name, this_chr):
    """
    Checks for exon, gene, and intron records from this chr
    """
    exon_select = """SELECT * from exon 
                     WHERE chr=\'%s\'
                     LIMIT 10""" % this_chr

    exon_records = db.getDBRecords_Dict(exon_select, db_name)

    if exon_records == []:
        return False
    
    intron_select = """SELECT * from intron 
                     WHERE chr=\'%s\'
                     LIMIT 10""" % this_chr

    intron_records = db.getDBRecords_Dict(intron_select, db_name)

    if intron_records == []:
        return False

    gene_select = """SELECT * from gene 
                     WHERE chr=\'%s\'
                     LIMIT 10""" % this_chr

    gene_records = db.getDBRecords_Dict(gene_select, db_name)

    if gene_records == []:
        return False

    # passed all checks
    return True

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getChr(tmp_file):
    tmp_file = tmp_file.rstrip(".gtf")

    return tmp_file.split("_")[-1]
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
