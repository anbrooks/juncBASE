# sequence_wrap.py
# Author: Angela Brooks
# Program Completion Date:
# Description: A wrapper class to reduce coding when trying to use dna sequences
# using Biopython
# Modification Date(s):
# 
# class DNA:
#    __init__(self, sequence)
#    __str__(self)
#    __add__(self, other)
#    __getslice__(self, i,j)
#    __getitem__(self,i)
#    __len__(self)
#    __cmp__(self, other)
#    tostring(self)
#    prettyString(self)
#   extractSeq(self,start,end)
#    reverseComplement(self)
#    find(self, other)
#    removeN()
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

import sys, getopt, pdb

from string import maketrans

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
###########
#CONSTANTS#
###########

DNA_ALPHABET = IUPAC.unambiguous_dna
    
###############
#END CONSTANTS#
###############


#########
#CLASSES#
#########
class DNA:
    def __init__ (self, sequence):
        self.seq = Seq( sequence, DNA_ALPHABET )
    def __str__ (self):
        return self.seq.tostring()
    def __add__ (self, other):
        if isinstance(other,DNA):
            return self.__class__(self.seq.__add__(other.seq).tostring())
        elif type(other) == type(""):
            return self.__class__(self.seq.__add__(DNA(other)).tostring())
        else:
            raise TypeError, ("adding incompatible types")
    def __getslice__(self, i,j):
        return self.__class__(self.seq[i:j].tostring())
    def __getitem__(self,i):
        return self.seq[i]
    def __len__(self):
        return len(self.seq)
    def __cmp__(self, other):
        if isinstance(other,DNA):
            return cmp(self.seq.tostring().upper(),other.seq.tostring().upper())
        elif type(other) == type(""):
            return cmp(self.seq.tostring().upper(),other.upper())
        else:
            raise TypeError, ("comparing incompatible types")

    
    def tostring (self):
        return self.seq.tostring()

    def prettyString(self, lineLen=50):
        """
        Many times one will want to print a sequence in a format with line
        breaks. This function will do so.
        """
        s = self.tostring()
        ps = ""
        
        for k in range(1,len(s) + 1):
            ps += s[k-1]
            if k % lineLen == 0:
                ps += "\n"

        # Make sure last char is \n
        if not ps.endswith("\n"):
            ps += "\n"

        return ps 

    # Used to extract a subsequence of the DNA sequence.  The start and end are
    # indexed starting with 1, like real sequences are.
    def extractSeq(self, start, end):
        return self.__class__(self.seq[start-1:end].tostring())

    def reverseComplement(self):
        l = list(self.seq.tostring().translate(maketrans("atcgATCG","tagcTAGC")))    
        l.reverse()
        return self.__class__(''.join(l))
    def find (self, other):
        if type(other) == type(""):
            return self.seq.tostring().find(other)

        return self.seq.tostring().find(other.seq.tostring())
    def removeN (self):
        seq_str = self.seq.tostring()
        seq_str = seq_str.replace("N","")
        seq_str = seq_str.replace("n","")
        return self.__class__(seq_str)    


class Locus:
    def __init__(self, seqname, start, end, source="Locus_object",
                 feature=None, score=None, strand=None, frame=None,
                 attributes=None):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame 
        self.attributes = attributes # This will likely be a gene name

    def __str__(self):
        """
        A Locus will basically print out in GFF format
        """
        strand = self.strand
        if self.strand is None:
            strand = "." 
    
        feature = self.feature
        if self.feature is None:
            feature = ""

        attributes = self.attributes
        if self.attributes is None:
            attributes = ""
 
        s = "%s\t%s\t%s\t%d\t%d\t" % (self.seqname,
                                      self.source,
                                      feature,
                                      self.start,
                                      self.end)
        if not self.score is None:
            s += "%f\t" % self.score
        else:
            s += "\t"

        s += "%s\t" % strand

        if not self.frame is None:
            s += "%d\t" % self.frame
        else:
            s += "\t"

        s += "%s" % self.attributes
                                     
        return s
#############
#END CLASSES#
#############
 
######
#MAIN#    
######
def main():
    seq1 = DNA("GATACA")
    print "The sequence is: %s" % seq1

    print "The sequences at positions 2-4 of",
    print "seq1 is: %s" % seq1.extractSeq(2,4)

    seq2 = DNA("TTTT")
    seq3 = seq1 + seq2
    print "Adding to the sequence: %s" % seq3

    seq4 = seq1 + DNA("GTATAT")
    print "Adding something different: %s" % seq4

    seq4 = seq4 + DNA("TTTC")
    print "Adding even more to previous: %s" % seq4

    print "length is: %d" % len(seq4)

    seq5 = seq4[0:3]
    print "The first codon is: %s" % seq5
    print type(seq5)

    seq6 = seq4.reverseComplement()
    print "The reverse complement is: %s" % seq6    

    str_sequence = seq6.tostring()
    print "Now it is a string: %s" % str_sequence
    print "See this is why: %s" % str_sequence.lower()    

    motif = DNA("TACT")
    print "The motif was found at: %d" % seq6.find(motif)

    seq7 = DNA("TTAGT")
    seq8 = DNA("TTAGt")

    print seq7 == seq8    

    find_str = "TACT"
    print "The string motif was found at:%d" % seq6.find(find_str)

##########
#END_MAIN#
##########

###########
#FUNCTIONS#
###########
def formatLine( line ):
    #format line
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

###############
#END FUNCTIONS#    
###############    
if __name__ == "__main__": main()
