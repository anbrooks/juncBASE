#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# runCuffmerge
# Author: Angela Brooks
# Program Completion Date:
# Description: Will clean up results from cufflinks denovo runs and combine into
# final gtf using cuffmerge
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb

import random

#############
# CONSTANTS #
#############
GFFREAD = "~/src/Cufflinks/cur/gffread"
CUFFMERGE = "~/src/Cufflinks/cur/cuffmerge"
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
    opt_parser.add_option("-f",
                          dest="file",
                          type="string",
                          help="""File containing path of files that will be
                                  merged into one gtf.""",
                          default=None)
    opt_parser.add_option("--no_single_exons",
                          dest="no_single_exons",
                          action="store_true",
                          help="""Will remove single exons that are de novo
                                  predictions from Cufflinks. Annotated single
                                  exons will remain.""",
                          default=False)
    opt_parser.add_option("-g",
                          dest="annotation",
                          type="string",
                          help="Annotation gtf to serve as a guide",
                          default=None)
    opt_parser.add_option("-o",
                          dest="cuffmerge_out",
                          type="string",
                          help="Directory of cuffmerge output. Def=\".\"",
                          default=".")
    opt_parser.add_option("-s",
                          dest="fasta_sequence",
                          type="string",
                          help="Location of genome sequence in fasta format.",
                          default=None)
    opt_parser.add_option("--tmp_dir",
                          dest="tmp_dir",
                          type="string",
                          help="""Temporary directory. Only required if using
                                  --no_single_exons""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-f")
    opt_parser.check_required("-g")
    opt_parser.check_required("-o")
    opt_parser.check_required("-s")

    gtf_files_file = options.file
    no_single_exons = options.no_single_exons
    if no_single_exons:
        opt_parser.check_required("--tmp_dir")
        tmp_dir = formatDir(options.tmp_dir)
        

    annotation = options.annotation
    cuffmerge_out = formatDir(options.cuffmerge_out)

    seq = options.fasta_sequence

    if no_single_exons:
        gtf_files = open(gtf_files_file)
        gtf_file_lines = gtf_files.readlines()
        gtf_files.close()

        gtf_files_for_cuffmerge = "%s/tmp_gtf_files.txt" % tmp_dir
        gtf_files_for_cuffmerge_fh = open(gtf_files_for_cuffmerge, "w")
                                      
        for line in gtf_file_lines:
            line = formatLine(line)
            # Random file name generated to prevent conflicts
            cleaned_file = "%s/cufflinks_%d_cleaned.gtf" % (tmp_dir,
                                                            random.randint(1,10000000000))
            cmd = "%s %s -o %s -T -U" % (GFFREAD,line, cleaned_file)
            os.system(cmd)
            gtf_files_for_cuffmerge_fh.write(cleaned_file + "\n")

        gtf_files_for_cuffmerge_fh.close()

    else:
        gtf_files_for_cuffmerge = gtf_files_file

    # Now run cuffmerge
    cmd = "%s -g %s " % (CUFFMERGE, annotation)
    cmd += "-o %s " % cuffmerge_out
    cmd += "-s %s " % seq
    cmd += gtf_files_for_cuffmerge
    os.system(cmd)
			
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

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
