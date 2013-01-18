#!/usr/bin/env python
# Co-written by Angadhjot Hundal.
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

import time
import re
import os
import sys
import fileinput
import pdb

KNOWN_TAGS = ["XA", "XM", "MD", "NM", "NH"]
JUNC_DELIM = "|"

def JUNCTION( w2 ):
    junc_read = w2.split("|")
    return junc_read

def READ( junc ):
    t = junc[-1]
    genes = junc[-3]
    chr = junc[0]
    s = junc[-2]
    intron_start = int(junc[2]) + 1
    intron_end = int(junc[3]) - 1

    return t, genes, chr, s, intron_start, intron_end

def INTRON( w2 ):
    temp = re.findall('\:\:\:\d+', w2)
    intron_end = temp[-1][3:]
    return intron_end

def STRAND_DIR( junc ):
        if junc == '-':
                direction = 16
        else:
                direction = 0
        return direction

def DIR( w1, dir ):

    if dir == 16:
        # Bitwise exclusive or when the junction is on the opposite
        # strand
        return repr(int(w1)^16)

    return w1

def CIGAR( w3 , w5, dir , intron_start , intron_end, sub ):
    total = int(w5[:-1])
    length = total - int(sub)
    intron_len = intron_end_pos - intron_start_pos + 1

    if dir == 0:
        first_M = length - int(w3) + 1
        sec_M = total - first_M
    else: 
        sec_M = length - int(w3) + 1
        first_M = total - sec_M
        
    exon_start_pos1 = int(intron_start) - first_M
    exon_end_pos1 = int(intron_start) - 1
    exon_start_pos2 = int(intron_end) + 1
    exon_end_pos2 = int(exon_start_pos2) + sec_M
    w3 = str(exon_start_pos1)
    w5 = str(first_M)+"M"+str(intron_len)+"N"+str(sec_M)+"M"

    return w3, w5

def SEQ( w9 ):
    reverse = []
    for base in w9:
        if base == "C":
            reverse.append("G")
        elif base == "G":
            reverse.append("C")
        elif base == "A":
            reverse.append("T")
        elif base == "T":
            reverse.append("A")
    reverse_compliment = reverse[::-1]
    final_seq = ''.join( reverse_compliment )
    return final_seq
        
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line


def reverseQual(old_qual):
    qual_list = list(old_qual)
    qual_list.reverse()
    return "".join(qual_list)

def reverseMD(value):
    md_list = re.split("(\d+)", value)
    md_list.reverse()
    return "".join(md_list)
    
def diffJunc(other_junc_format):

        (chr, upstr_start, upstr_end, downstr_start, downstr_end, type, strand, 
         gene1, gene2) = other_junc_format.split(":::")

        # THIS IS FOR A REALLY WEIRD FORMAT
        intron_start = int(upstr_end) - 1

        new_format = [type, gene1, gene2, chr.lstrip("chr"),strand,
                      repr(intron_start), downstr_start]
        return ":::".join(new_format)

usage = "usage: %s infile outfile overhang [--watsonJunctions]" % os.path.basename(sys.argv[0])

if len( sys.argv ) < 4:
    print usage
    sys.exit(1)
else:
#    print "There are %s args " %len( sys.argv )

    fileToSearch      = open( sys.argv[1], 'r' )
    fileToOutput      = sys.argv[2]
    oldFileName  = 'old_' + fileToOutput
    tempFileName = 'temp_' + fileToOutput
# If there's already an 'old_' prefixed backup file there from
# a previous run, remove it...
    if os.path.isfile( tempFileName ):
        os.remove( tempFileName )
    if os.path.isfile( oldFileName ):
        os.remove( oldFileName )
    if os.path.isfile( fileToOutput ):
        os.remove( fileToOutput )    

    sub = sys.argv[3]

    # Junction sequences can be given in the watson strand only, which will
    # affect the way the junction is converted.
    isWatsonJunctions = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "--watsonJunctions":
            isWatsonJunctions = True
        

    # Remove unaligned reads
    fTemp = open(tempFileName , 'a')
#    fTemp.writelines("This is a temporary file: \n")    
    for lines in fileToSearch:
        words = lines.split()
        if words[1] != '4':
            fTemp.write(lines)

    fileToSearch.close()
    fTemp.close()

    fin = open(tempFileName , 'r')
    fout = open(oldFileName , 'a')

    insertionFlag = False
    deletionFlag = False
    softFlag = False
    hardFlag = False
    
    for lines in fin:
        lines = formatLine(lines)
        words = lines.split()
        if "I" in words[5]:
            if not insertionFlag:
                print "Not handling insertion alignments yet..."
                insertionFlag = True
            continue
        if "D" in words[5]:
            if not deletionFlag:
                print "Not handling deletion alignments yet..."
                deletionFlag = True
            continue
        if "S" in words[5]:
            if not softFlag:
                print "Not handling soft clipping alignments yet..."
                softFlag = True
            continue
        if "H" in words[5]:
            if not hardFlag:
                print "Not handling hard clipping alignments yet..."
                hardFlag = True
            continue

        if len(words) >= 11:
####        if words[2].startswith("chr"):
####            words[2] = diffJunc(words[2])
            junc_read = JUNCTION( words[2] )
            if len(junc_read) > 1: # It is a junction read.
#                type, gene1, gene2, read, strand, intron_start_pos = READ( junc_read )
                type, genes, chr, strand, intron_start_pos, intron_end_pos = READ( junc_read )
                
#                direction = STRAND_DIR( strand )
                direction = int(words[1])
                new14 = "Y0:Z:"+words[2]
                new15 = "XS:A:"+strand
                if not isWatsonJunctions:
                    words[1] = DIR( words[1], direction )
                    if strand == "-" and direction == 0:
                        words[9] = SEQ( words[9] )
                        words[10] = reverseQual(words[10])
                words[2] =  chr
                words[3], words[5] = CIGAR( words[3] , words[5], int(words[1]) , intron_start_pos , intron_end_pos , sub)
                remaining_fields = words[11:]
                new_fields = []
                for field in remaining_fields:
                    tag, type, val = field.split(":")
                    if tag not in KNOWN_TAGS:
                        print "ERROR: Unknown tag: %s" % tag
                        sys.exit(1)
                    if tag == "XA" or tag == "XM" or tag == "NM" or tag == "NH":
                        new_fields.append(field)
                    if tag == "MD":
                        if strand == "-":
                            new_field = "%s:%s:%s" % (tag,
                                          type,
                                          reverseMD(val))
                            new_fields.append(field)
                        else:
                            new_fields.append(field)
                
                    
#                t1 = words[11]
#                        words[11] = words[13]
#                            words[13] = "X"+t1[-1:]+":i:1"
                new_fields.append(new14)
                new_fields.append(new15)
                new_line = [ words[0], words[1], words[2], words[3], words[4], words[5], words[6], words[7], 
                         words[8], words[9], words[10]]
                new_line.extend(new_fields)
                new_string = '\t'.join(new_line)
                fout.write(new_string+'\n')                
            else: # Is a genome read.
                fout.write(lines + "\n")
                break
        else:
            print "SAM file needs at least 11 fields."
            sys.exit(1)

    fin.close()
    fout.close()

# Rename the original file by prefixing it with 'old_'
os.rename( oldFileName, fileToOutput )

# Remove temporary file
if os.path.isfile( tempFileName ):
    os.remove( tempFileName )

