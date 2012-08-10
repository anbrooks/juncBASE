#!/lab/64/bin/python
# createPseudoSample.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.
"""Creates input files for a pseudo sample. Used when determining alternatively
spliced events with more than two samples.
"""

import sys
import optparse 
import os
import gzip
import pdb

from helperFunctions import runCmd
#############
# CONSTANTS #
#############
SHELL = "/bin/tcsh"

NUM_GENOME_READS = 10
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
                          dest="input_dir",
                          type="string",
                          help="""Root directory containing all input
                                  directories that will be used for
                                  getASEventReadCounts.py""",
                          default=None)
    opt_parser.add_option("-s",
                          dest="sample",
                          type="string",
                          help="""Real sample from which the pseudo sample will be
                                  derived from. This is an arbitrary decision
                                  and the choice should have no effect on the
                                  eventual result.""",
                          default=None)
    opt_parser.add_option("--by_chr",
                          dest="by_chr",
                          action="store_true",
                          help="""Will indicate that the input files were split by
                                  chr""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-s")

    by_chr = options.by_chr

    input_dir = None
    if os.path.exists(options.input_dir):
        input_dir = os.path.abspath(options.input_dir)
    else:
        print "Input directory does not exist."
        opt_parser.print_help()
        sys.exit(1)
    if input_dir.endswith("/"):
        input_dir = input_dir.rstrip("/")

    samp = options.sample

    # if by_chr is False, only one chr
    chr2pseudo_dir = {}
    
    # Create pseudo directory
    if by_chr:
        # grab chr names from preProcess.bed files
        chrs = getChr(input_dir)
        if chrs == []:
            print "Missing bed files from preProcess_getASEventReadCounts_by_chr_step2.py. They should be in: %s" % input_dir
            opt_parser.print_help()
            sys.exit(1)

        for chr in chrs:
            chr2pseudo_dir[chr] = input_dir + "/pseudo_" + chr
    else:
        pseudo_dir = input_dir + "/pseudo"
        chr2pseudo_dir["all"] = pseudo_dir

    for chr in chr2pseudo_dir:
        if os.path.exists(chr2pseudo_dir[chr]):
            # Remove existing files
            for subfile in os.listdir(chr2pseudo_dir[chr]):
                try:
                    os.remove(chr2pseudo_dir[chr] + "/" + subfile)
                except:
                    pdb.set_trace()
        else:
            os.mkdir(chr2pseudo_dir[chr])

    # Make a link to the junction file made in
    # preProcess_getASEventReadCounts_step2.py
    if "all" in chr2pseudo_dir:
        bed_file = input_dir + "/tmp_preProcess_getASEventReadCounts_step2.bed"
        if not os.path.exists(bed_file):
            print "Missing %s. Please run preProcess_getASEventReadCounts_step2.py" % bed_file
            sys.exit(1)

        ln_cmd = "ln -s %s %s/pseudo_junctions.bed" % (bed_file,
                                                       chr2pseudo_dir["all"])
        os.system(ln_cmd)

    else:
        for chr in chr2pseudo_dir:
            bed_file = input_dir + "/tmp_" + chr + "_preProcess_getASEventReadCounts_step2.bed"
            ln_cmd = "ln -s %s %s/pseudo_%s/pseudo_%s_junctions.bed" % (bed_file,
                                                                        input_dir,
                                                                        chr,
                                                                        chr)
            os.system(ln_cmd)

    # Link to intron_exon_junction file
    if "all" in chr2pseudo_dir:
        i_e_file = input_dir + "/%s/%s_intron_exon_junction_counts.txt" % (samp,
                                                                           samp)
        if not os.path.exists(i_e_file):
            print "Missing file %s." % i_e_file
            print "Please run preProcess_getASEventReadCounts_step3.py on sample %s." % samp
            sys.exit(1)

        ln_cmd = "ln -s %s %s/pseudo_intron_exon_junction_counts.txt" % (i_e_file,
                                                                         chr2pseudo_dir["all"])
        os.system(ln_cmd)

    else:
        for chr in chr2pseudo_dir:
            i_e_file = input_dir + "/%s/%s_%s/%s_%s_intron_exon_junction_counts.txt" % (samp,
                                                                                        samp, chr,
                                                                                        samp, chr)
            if not os.path.exists(i_e_file):
                print "Missing file %s." % i_e_file
                print "Please run preProcess_getASEventReadCounts_step3.py on sample %s and chr %s." % (samp,
                                                                                                        chr)
                sys.exit(1)

            ln_cmd = "ln -s %s %s/pseudo_%s_intron_exon_junction_counts.txt" % (i_e_file,
                                                                                chr2pseudo_dir[chr],
                                                                                chr)
            os.system(ln_cmd)
                    

    # Copy over a little of the genome reads
    if "all" in chr2pseudo_dir:
        genome_file = input_dir + "/%s/%s_genome_reads.txt.gz" % (samp, samp)
        pseudo_genome_file = pseudo_dir + "/pseudo_genome_reads.txt.gz"
        
        if not os.path.exists(genome_file):
            print "Missing file %s" % genome_file
            sys.exit(1)

        genome_fh = gzip.open(genome_file, "rb")
        pseudo_fh = gzip.open(pseudo_genome_file, "wb")

        line_ctr = 0
        for line in genome_fh:
            if NUM_GENOME_READS == line_ctr:
                break

            pseudo_fh.write(line)

            line_ctr += 1
        
        genome_fh.close()
        pseudo_fh.close()

    else:                                                                        
        for chr in chr2pseudo_dir:
            genome_file = input_dir + "/%s/%s_%s/%s_%s_genome_reads.txt.gz" % (samp,
                                                                              samp, chr,
                                                                              samp, chr)
            pseudo_file = chr2pseudo_dir[chr] + "/pseudo_%s_genome_reads.txt.gz" % chr

            if not os.path.exists(genome_file):
                print "Missing file %s" % genome_file
                sys.exit(1)

            genome_fh = gzip.open(genome_file, "rb")
            pseudo_fh = gzip.open(pseudo_file, "wb")

            line_ctr = 0
            for line in genome_fh:
                if NUM_GENOME_READS == line_ctr:
                    break
   
                pseudo_fh.write(line)
                line_ctr += 1

            genome_fh.close()
            pseudo_fh.close() 

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

def getChr(input_dir):
    chr_list = []
    for subfile in os.listdir(input_dir):
        if not subfile.endswith("preProcess_getASEventReadCounts_step2.bed"):
            continue

        first_split = subfile.split("preProcess")[0]
        second_split = first_split.split("tmp")[-1]
    
        chr = second_split.strip("_")

        if chr == "":
            continue

        chr_list.append(second_split.strip("_"))

    return chr_list
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
