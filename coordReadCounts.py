#!/usr/bin/env python
# coordReadCounts.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.
"""Outputs read counts for specified start and end coordinates.
    Uses AVL trees for faster searching.
"""
import sys
import optparse 
import re
import pdb
import gzip


from opus7.avlTree import AVLTree
############
# CONSTANTS #
#############
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
    opt_parser.add_option("--reads",
                          dest="read_file",
                          type="string",
                          help="Contains coordinates of reads. query_name, chr, start, end",
                          default=None)
    opt_parser.add_option("--sam",
                          dest="sam_file",
                          type="string",
                          help="""Contains coordinates of read alignments in sam
                                  format.""",
                          default=None)
    opt_parser.add_option("-l",
                          dest="read_len",
                          type="int",
                          help="""Length of the reads. Option given only if reads
                                  are of the same length.""",
                          default=None)
    opt_parser.add_option("--coords",
                          dest="coord_file",
                          type="string",
                          help="Optional: Columns are chr, start, end",
                          default=None)
    opt_parser.add_option("--sqlite_db_dir",
                          dest="sqlite_db_dir",
                          type="string",
                          help="""Location of sqlite database. If sqlite database
                                  is used, will override usage of MySQL database.""",
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
    opt_parser.add_option("--db_name",
                          dest="db_name",
                          type="string",
                          help="Optional: Get coordinates from this database.",
                          default=None)
    opt_parser.add_option("--db_table",
                          dest="table",
                          type="string",
                          help="""Optional: Get coordinates from this table,
                                  from specified database.  The table should
                                  have a chr, start, end""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output",
                          type="string",
                          help="Ouptut file of chr, start, end, and counts",
                          default=None)
    opt_parser.add_option("--read_assoc",
                          dest="read_assoc_file",
                          type="string",
                          help="""Ouptut of read coordinate associated with the
                                coord coordinate""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--reads")
    opt_parser.check_required("-o")

    if (not options.read_file) and (not options.sam_file):
        print "Must provide a read file."
        sys.exit(1)

    isSam = False
    if options.read_file:
        if options.read_file.endswith(".gz"):
            read_file = gzip.open(options.read_file, "rb")
        else:
            read_file = open(options.read_file)
    else:
        read_file = open(options.sam_file)
        isSam = True

    output = open(options.output, "w")

    read_assoc_file = None
    if options.read_assoc_file:
        read_assoc_file = open(options.read_assoc_file, "w")

    read_len = options.read_len

    if options.coord_file is None and options.db_name is None:
        print "Need to specify a coordinate file or database name and table."
        opt_parser.print_help()
        sys.exit(1)

    if options.coord_file and options.db_name:
        print "Only specify a coordinate file OR a database name and table."
        opt_parser.print_help()
        sys.exit(1)

    if options.db_name and options.table is None:
        print "Need to specify a table name."
        opt_parser.print_help()
        sys.exit(1)
   
    # Status
    print "Creating tree"
    if options.coord_file:
        coord_file = open(options.coord_file)

        chr2coord2count, chr2coordTree = parseCoordFile(coord_file)

    else:
        db = None
        if options.sqlite_db_dir:
            import pysqlite_wrap
            db = pysqlite_wrap.DB(options.sqlite_db_dir)
        else: # Use MySQL db
            import mysqldb_wrap 
            db = mysqldb_wrap.DB(options.host, options.user, options.passwd)   
     
        db_name = options.db_name
        table = options.table

        chr2coord2count, chr2coordTree = getCoordFromDB(db, db_name, table)


    # Status
    print "Processing Reads"
    for line in read_file:
        line = formatLine(line)

        if line.startswith("#"):
            continue
        if line.startswith("@"):
            continue

        if isSam:
            q_name, chr, start, end = parseSam(line)

            if q_name is None:
                continue
        else:
            line_list = line.split("\t")
            q_name = line_list[0]
            chr = line_list[1]
            start = int(line_list[2])

            try:
                end = int(line_list[3])
            except:
                if read_len:
                    end = start + read_len - 1
                else:
                    print "Expecting read length: %d" % read_len
                    sys.exit(1)

        if chr.startswith("chr"):
            chr = chr.replace("chr","")    

        if chr not in chr2coord2count:
            continue

        read_coord = (start, end)

        # Recursive function.  The tree is used for faster searching, the other
        # dictionary is used to store counts.
        if read_assoc_file:
            findReadHits(chr, q_name, read_coord, chr2coord2count[chr],
                         chr2coordTree[chr], read_assoc_file)
        else:
            findReadHits(chr, q_name, read_coord, chr2coord2count[chr], chr2coordTree[chr])

    # Print to output
    for chr in chr2coord2count:
        for coord in chr2coord2count[chr]:
            out_str = "%s\t%d\t%d\t%d\n" % (chr,
                                            coord[0],
                                            coord[1],
                                            chr2coord2count[chr][coord])
            output.write(out_str)
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def findReadHits(chr, q_name, read_coord, coord2count, coordTree, read_assoc_file=None):
    """
    Uses a binary tree to traverse through input coordinates in order to find
    reads that are completely contained into the coordinate. Increments counts
    for the coord in the other dictionary
    """
    this_tree = coordTree
    thisIsContained = False
    if isContained(read_coord, this_tree.getKey()):
        # Add counts
        coord2count[this_tree.getKey()] += 1
    
        # Print to association file if necessary
        if read_assoc_file:
            out_line = "%s\t%s:%d-%d\t%s:%d-%d\n" % (q_name,
                                                     chr,
                                                     read_coord[0],
                                                     read_coord[1],
                                                     chr,
                                                     this_tree.getKey()[0],
                                                     this_tree.getKey()[1])
            read_assoc_file.write(out_line)

        thisIsContained = True

    if not this_tree.getIsLeaf():
        if read_coord < this_tree.getKey() or thisIsContained:
            # find Hits in all subTrees if read_coord checks out
            if moveLeftCheck(read_coord, this_tree.getLeft()):
                findReadHits(chr, q_name, read_coord, coord2count, this_tree.getLeft(),
                             read_assoc_file)

        if read_coord > this_tree.getKey() or thisIsContained:
            # find Hits in all subTrees if read_coord checks out
            if moveRightCheck(read_coord, this_tree.getRight()):
                findReadHits(chr, q_name, read_coord, coord2count,
                             this_tree.getRight(), read_assoc_file)

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def formatChr(chr):
    if chr.startswith("chr"):
        return chr.replace("chr","")

    return chr

def getCoordFromDB(db, db_name, table):
    """ 
    Returns two dictionaries:
    chr2coord2count: Holds the counts for each coordinate
    chr2coordTree: Holds each search tree for each chr in the coord file
    """
    chr2coord2count = {}
    chr2coordTree = {}

    select_state = "SELECT chr, start, end FROM %s" % table
    records = db.getDBRecords_Dict(select_state, db_name)

    for row in records:
        chr = row["chr"]
        start = int(row["start"])
        end = int(row["end"])

        coord = (start, end)

        if chr.startswith("chr"):
            chr = chr.replace("chr","")

        if chr in chr2coord2count:
            chr2coord2count[chr][coord] = 0
        else:
            chr2coord2count[chr] = {coord:0}

        if chr in chr2coordTree:
            chr2coordTree[chr].insert(coord)
        else:
            chr2coordTree[chr] = AVLTree()
            chr2coordTree[chr].insert(coord)

    return chr2coord2count, chr2coordTree

def isContained(first_coord, second_coord):
    if first_coord == second_coord:
        return True

    if first_coord[0] >= second_coord[0]:
        if first_coord[1] <= second_coord[1]:
            return True

    return False

def moveLeftCheck(read_coord, this_tree):
    """
    If moving left, the read_coord should be contained in or less than the left coord.
    If these don't hold, then the left_tree has to have right children
    """
    if this_tree.isEmpty:
        return False

    if read_coord < this_tree.getKey():
        return True

    if isContained(read_coord, this_tree.getKey()):
        return True

    # Check for right children
    if not this_tree.getRight().isEmpty:
        return True

    return False

def moveRightCheck(read_coord, this_tree):
    """
    If moving right, the read_coord should be contained in or more than the
    right coord.
    If these don't hold, then the right_tree has to have left children
    """
    if this_tree.isEmpty:
        return False

    if read_coord > this_tree.getKey():
        return True

    if isContained(read_coord, this_tree.getKey()):
        return True

    # Check for left children
    if not this_tree.getLeft().isEmpty:
        return True
    
    return False

def parseCoordFile(coord_file):
    """ 
    Returns two dictionaries:
    chr2coord2count: Holds the counts for each coordinate
    chr2coordTree: Holds each search tree for each chr in the coord file
    """
    chr2coord2count = {}
    chr2coordTree = {}

    for line in coord_file:
        line = formatLine(line)

        line_list = line.split("\t")
        chr = line_list[0] 
        start = int(line_list[1])
        end = int(line_list[2])

        coord = (start, end)

        if chr.startswith("chr"):
            chr = chr.replace("chr","")

        if chr in chr2coord2count:
            chr2coord2count[chr][coord] = 0
        else:
            chr2coord2count[chr] = {coord:0}

        if chr in chr2coordTree:
            chr2coordTree[chr].insert(coord)
        else:
            chr2coordTree[chr] = AVLTree()
            chr2coordTree[chr].insert(coord)

    return chr2coord2count, chr2coordTree

def parseSam(line, maxGap=20):
    """
    Returns the following:
    q_name, chr, start, end

    Will ignore reads with gaps larger than the maxGap.
    """
    line_list = line.split("\t")
    q_name = line_list[0]

    chr = line_list[2]
    chr = formatChr(chr)    

    start = int(line_list[3])

    cigar = line_list[5]

    # Will not deal with these types of alignments for now
    if ("N" in cigar) or ("S" in cigar) or ("H" in cigar) or ("P" in cigar):
        return None, None, None, None

    total_len = 0
    matches = re.findall("\d+M", cigar)
    for m in mathces:
        this_len = int(m.rstrip("M"))
        total_len += this_len
    insertions = re.findall("\d+I", cigar)
    for i in matches:
        this_len = int(i.rstrip("I"))
        if this_len > maxGap:
            return None, None, None, None
        total_len += this_len
    
       
    return q_name, chr, start, start + total_len - 1
    
        

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
