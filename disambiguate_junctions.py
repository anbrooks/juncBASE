#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# disambiguate_junctions.py
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
import shutil

from Bio import SeqIO
from getASEventReadCounts import convertCoordStr
#from createPseudoSample import getChr
#############
# CONSTANTS #
#############

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
                          help="Directory of juncBASE input files",
                          default=None)
    opt_parser.add_option("-g",
                          dest="genome_fasta",
                          type="string",
                          help="""Fasta file of genome sequence. Can also be one
                                  or a few chromosomes.""",
                          default=None)
    opt_parser.add_option("--by_chr",
                          dest="by_chr",
                          action="store_true",
                          help="""Indicates that input files are organized by
                                  chromosome""",
                          default=False)
    opt_parser.add_option("--majority_rules",
                          dest="majority_rules",
                          action="store_true",
                          help="""Will ensure that strand of each junction across all samples
                                  are the same using the strand that is most
                                  aligned to. This option should be used if
                                  there are junctions corresponding to
                                  non-major and non-minor spliceosome splice site sequences.""",
                          default=False)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-g")

    majority_rules = options.majority_rules

    input_dir = formatDir(options.input_dir)

    samples = getSamples(input_dir)

    genome_file = options.genome_fasta
    
    by_chr = options.by_chr

    try:
        records = SeqIO.index(genome_file, "fasta")
    except:
        print "Could not open genome sequence."
        sys.exit(1)

    if by_chr:
        chr_list = getChr(input_dir, samples)

        for this_chr in chr_list:
            if this_chr in records:
                chr_seq = records[this_chr].seq
            elif unFormatChr(this_chr) in records:
                chr_seq = records[unFormatChr(this_chr)].seq
            else:
                print "Cannot find %s in genome sequence." % this_chr
                sys.exit(1)
            jcn2strands = {} 
            for this_samp in samples:
                bed_file_name = "%s/%s/%s_%s/%s_%s_junctions.bed" % (input_dir,
                                                                     this_samp,
                                                                     this_samp, this_chr,
                                                                     this_samp, this_chr)

                # Copy old bed file to keep raw information.
                original_bed_name = "%s/%s/%s_%s/%s_%s_junctions_original.bed" % (input_dir,
                                                                                  this_samp,
                                                                                  this_samp, this_chr,
                                                                                  this_samp, this_chr)
                # Do not overwrite original file if it exists (prevents totally
                # losing information in reruns
                if not os.path.exists(original_bed_name):
                    shutil.copy(bed_file_name, original_bed_name)
                
                original_bed = open(original_bed_name)
                bed_file = open(bed_file_name, "w")
    
                for line in original_bed:
                    if line.startswith("track"):
                        bed_file.write(line)
                        continue

                    line_list = line.split("\t")
    
                    if line_list[5] == ".":
                        # Updates the strand position in the list
                        disambiguateJcnStr(chr_seq, line_list, majority_rules) 

                    if majority_rules:
                        updateJcn2Strands(jcn2strands, 
                                          line_list[3], 
                                          line_list[5])
                    
                    bed_file.write("\t".join(line_list)) 

                original_bed.close()
                bed_file.close()

            # Fix junctions across all samples
            if majority_rules:
                fixed_jcn2strand = getMajorityStrand(jcn2strands)
    
                if fixed_jcn2strand != {}:
                    for this_samp in samples:
                        bed_file_name = "%s/%s/%s_%s/%s_%s_junctions.bed" % (input_dir,
                                                                             this_samp,
                                                                             this_samp, this_chr,
                                                                             this_samp, this_chr)

                        bed_file = open(bed_file_name)
                        bedlines = bed_file.readlines()
                        bed_file.close()

                        bed_file = open(bed_file_name, "w")
                        for line in bedlines:
                            if line.startswith("track"):
                                bed_file.write(line)
                                continue

                            line_list = line.split("\t")
                           
                            try: 
                                if line_list[3] in fixed_jcn2strand:
                                    line_list[5] = fixed_jcn2strand[line_list[3]]
                            except:
                                pdb.set_trace()
        
                            bed_file.write("\t".join(line_list))

                        bed_file.close()
                    

    # Files are just organized by sample
    else:
        jcn2strands = {}
        for this_samp in samples:
            bed_file_name = "%s/%s/%s_junctions.bed" % (input_dir,
                                                        this_samp,
                                                        this_samp)

            original_bed_name = "%s/%s/%s_junctions_original.bed" % (input_dir,
                                                                     this_samp,
                                                                     this_samp)

            # Do not overwirte original file
            if not os.path.exists(original_bed_name):
                shutil.copy(bed_file_name, original_bed_name)

            original_bed = open(original_bed_name)
            bed_file = open(bed_file_name, "w")

            for line in original_bed:
                if line.startswith("track"):
                    bed_file.write(line)
                    continue

                line_list = line.split("\t")

                if line_list[5] == ".":
                    disambiguateJcnStr_findChr(records, line_list, majority_rules)

                if majority_rules:
                    updateJcn2Strands(jcn2strands, 
                                      line_list[3], 
                                      line_list[5])

                bed_file.write("\t".join(line_list))

            original_bed.close()
            bed_file.close

        # Fix junctions across all samples
        if majority_rules:
            fixed_jcn2strand = getMajorityStrand(jcn2strands)

            if fixed_jcn2strand != {}:
                for this_samp in samples:
                    bed_file_name = "%s/%s/%s_junctions.bed" % (input_dir,
                                                        this_samp,
                                                        this_samp)
                    bed_file = open(bed_file_name)
                    bedlines = bed_file.readlines()
                    bed_file.close()

                    bed_file = open(bed_file_name, "w")
                    for line in bedlines:
                        if line.startswith("track"):
                            bed_file.write(line)
                            continue

                        line_list = line.split("\t")
                        
                        if line_list[3] in fixed_jcn2strand:
                            line_list[5] = fixed_jcn2strand[line_list[3]]
    
                        bed_file.write("\t".join(line_list))
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def disambiguateJcnStr(chr_seq, line_list, majority_rules):
    """
    Will use splice site sequence to infer strand
    """

    try:
        chr, start, end = convertCoordStr(line_list[3])
    except:
        print "Junction BED file must have intron position in 4th column."
        sys.exit(1)    

    intron_seq = chr_seq[start-1:end]

    if intron_seq.startswith("GT") and intron_seq.endswith("AG"):
        line_list[5] = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("AC"):
        line_list[5] = "-"
    # Other common splice site sequence
    elif intron_seq.startswith("GC") and intron_seq.endswith("AG"):
        line_list[5] = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("GC"):
        line_list[5] = "-"
    # minor spliceosome
    elif intron_seq.startswith("AT") and intron_seq.endswith("AC"):
        line_list[5] = "+"
    elif intron_seq.startswith("GT") and intron_seq.endswith("AT"):
        line_list[5] = "-"
    # Priority to 5' splice site since there is more information
    # there
    elif intron_seq.startswith("GT"):
        line_list[5] = "+"
    elif intron_seq.endswith("AC"):
        line_list[5] = "-"
    elif intron_seq.endswith("AG"):
        line_list[5] = "+"
    elif intron_seq.startswith("CT"):
        line_list[5] = "-"
    else:
        if not majority_rules: # Strand will resolved later if majority_rules
            print "Cannot find strand for %s" % line_list[3]


def disambiguateJcnStr_findChr(records, line_list, majority_rules):
    """
    Will identify chr sequence, then call disambiguateJcnStr
    """
    this_chr = line_list[0]
 
    chr_seq = None
    if this_chr in records:
        chr_seq = records[this_chr].seq
    elif unFormatChr(this_chr) in records:
        chr_seq = records[unFormatChr(this_chr)].seq
    else:
        print "Cannot find %s in genome sequence." % this_chr
        sys.exit(1)

    disambiguateJcnStr(chr_seq, line_list, majority_rules)
    

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getChr(input_dir, samples):
    """
    Returns a list of all chromosomes
    """
    chr_names = set([])

    for samp in samples:
        samp_dir = input_dir + "/" + samp
        for subdir in os.listdir(samp_dir):
            chr = subdir.replace(samp, "")
            chr = chr.lstrip("_") 
            chr_names.add(chr)

    return list(chr_names)

def getMajorityStrand(jcn2strands):
    fixed_jcn2strand = {}
    for jcn in jcn2strands:
        if len(jcn2strands[jcn]) == 1:
            if "." in jcn2strands[jcn]:
                print "No strand for %s" % jcn
            continue

        if "+" in jcn2strands[jcn]:
            pos_count = jcn2strands[jcn]["+"]
        else:
            pos_count = 0
        if "-" in jcn2strands[jcn]:
            neg_count = jcn2strands[jcn]["-"]
        else:
            neg_count = 0

        if pos_count == neg_count:
            fixed_jcn2strand[jcn] = "."
            print "No strand for %s" % jcn
            continue

        if pos_count > neg_count:
            fixed_jcn2strand[jcn] = "+"
        else:
            fixed_jcn2strand[jcn] = "-"

    return fixed_jcn2strand

def getSamples(input_dir):
    samples = []

    for this_samp in os.listdir(input_dir):
        if "pseudo" in this_samp:
            continue

        samp = input_dir + "/" + this_samp

        if not os.path.isdir(samp):
            continue

        samples.append(this_samp)

    return samples

def unFormatChr(chr_name):
    return chr_name.lstrip("chr")

def updateJcn2Strands(jcn2strands, jcn, strand):
    if jcn in jcn2strands:
        if strand in jcn2strands[jcn]:
            jcn2strands[jcn][strand] += 1
        else: 
            jcn2strands[jcn][strand] = 1
    else:
        jcn2strands[jcn] = {strand:1}
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
