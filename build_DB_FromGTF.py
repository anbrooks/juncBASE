#!/lab/64/bin/python
# build_DB_FromGTF.py
# Author: Angela Brooks
# Program Completion Date:
# Description: Takes .gtf files as input and parses out the information to
# build a database containing the transcript annotations 
# Modification Date(s):
# I highly recommend obtaining GTF files from ensembl.  It is important to have
# transcript ids associated with genes.
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb

from pysqlite_wrap import DB
from sequence_wrap import DNA

from Bio import SeqIO

from helperFunctions import updateDictOfLists

#############
# CONSTANTS #
#############
#VERSION = ""
# Database variables
DB_DIR = "./"

# TEST
#VALID_CHR = ["T"]
#VALID_CHR = ["4"]


#VALID_CHR = ["2L", "2LHet", "2R", "2RHet", "3L", "3LHet",
#             "3R", "3RHet", "4", "X", "XHet", "YHet", "U", "M"]

NOT_FOUND = -1
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

class GFF_Record:
    def __init__(self, tab_delim_line):
        #Error check
        if tab_delim_line.find("\t") == NOT_FOUND:
            print "GFF file must be tab delimited."
            sys.exit(1)

        l = tab_delim_line.split("\t")

        #Error check
        if len(l) > 9:
            print "Too many fields in GFF line:"
            print tab_delim_line
            sys.exit(1)

        #initialize variables
        self.seqname = ""
        self.source = ""
        self.feature = ""
        self.start = ""
        self.end = ""
        self.score = ""
        self.strand = ""
        self.frame = ""
        self.attributes = ""

        for i in range(len(l)):
            if i == 0:
                self.seqname = l[0]
            if i == 1:
                self.source = l[1]
            if i == 2:
                self.feature = l[2]
            if i == 3:
                self.start = l[3]
            if i == 4:
                self.end = l[4]
            if i == 5:
                self.score = l[5]
            if i == 6:
                self.strand = l[6]
                if ((self.strand != "+") and 
                    (self.strand != "-")):
                    self.strand = "."
                
            if i == 7:
                self.frame = l[7]
                if ((self.frame != 0) and (self.frame != 1) and 
                    (self.frame != 2)):
                    self.frame = "."
            if i == 8:
                self.attributes = l[8]

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
                          help="""Will initialize the database by creating the
                                  database and tables. Will not insert any
                                  data. If database already exists, it will
                                  erase the existing data.""",
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
    opt_parser.add_option("-f",
                          dest="genome_file_name",
                          type="string",
                          help="""Fasta file containing all chromosome
                                  sequences.  If this option is given, exon and
                                  intron sequences will be stored in the
                                  database as well. If running this script on
                                  individual chromosomes (for parallelization),
                                  only use give the chromosome sequence, not the
                                  full genome""",
                          default=None)
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
    

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-d")

    db_obj = DB(options.sqlite_db_dir)

    db_name = options.db_name

    initialize = options.initialize

    if initialize:
        db_obj.createDatabase(db_name)
        createAnnotTables(db_obj, db_name)
        sys.exit(0)

    # validate the command line arguments
    opt_parser.check_required("-g")

    use_gene_name = options.use_gene_name

    genome_file = options.genome_file_name

    chr2gtf_lines = getGTFLines(options.gtf_file)

    buildAnnotDB(db_obj, chr2gtf_lines, db_name, use_gene_name)
   
    if genome_file:
        extractInsertSeqs("exon", db_obj, genome_file, db_name)
        extractInsertSeqs("intron", db_obj, genome_file, db_name)
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildAnnotDB(db, chr2gtf_lines, db_name, use_gene_name):

    for chr in chr2gtf_lines:
        # Dictionary to hold information before creating exon entries.
        # Dictionaries will be of the form {transcript_id: [gff_obj]}
        five_utr_dict = {}
        three_utr_dict = {}
        start_codon_dict = {}
        stop_codon_dict = {}
        cds_dict = {}
        exon_dict = {}

        # Holds transcript and strand information separately
        # {transcript_id: strand}
        transcript_id2strand = {}

        for line in chr2gtf_lines[chr]:
            gff_obj = GFF_Record(line)

            if gff_obj.feature != 'exon':
                continue

            transcript_id = get_transcript_id(gff_obj) 
        
            # Update transcript_id2strand dictionary
            if not transcript_id in transcript_id2strand:
                transcript_id2strand[transcript_id] = gff_obj.strand

            # Update all feature dictionaries
#           if gff_obj.feature == '5UTR':
#               updateDictOfLists(five_utr_dict, transcript_id, gff_obj)
#           elif gff_obj.feature == '3UTR':     
#               updateDictOfLists(three_utr_dict, transcript_id, gff_obj)
#           elif gff_obj.feature == 'start_codon':     
#               updateDictOfLists(start_codon_dict, transcript_id, gff_obj)
#           elif gff_obj.feature == 'stop_codon':     
#               updateDictOfLists(stop_codon_dict, transcript_id, gff_obj)
#           elif gff_obj.feature == 'CDS':     
#               updateDictOfLists(cds_dict, transcript_id, gff_obj)
            if gff_obj.feature == 'exon':
                updateDictOfLists(exon_dict, transcript_id, gff_obj)
            else:
                print "Not using feature %s" % (gff_obj.feature)


        buildExonTable(db, chr, transcript_id2strand,
                       five_utr_dict,
                       three_utr_dict,
                       start_codon_dict,
                       stop_codon_dict,
                       cds_dict,
                       exon_dict, db_name, use_gene_name)

#       insertFeature(db, "cds", cds_dict, db_name)
#       insertFeature(db, "start_codon", start_codon_dict, db_name)
#       insertFeature(db, "stop_codon", stop_codon_dict, db_name)
#       insertFeature(db, "five_utr", five_utr_dict, db_name)
#       insertFeature(db, "three_utr", three_utr_dict, db_name)

        buildGeneTable(db, db_name, chr)       

        inferIntrons(db, db_name, chr)
 
def buildExonTable(db, chr, transcript_id2strand,
                   five_utr_dict,
                   three_utr_dict,
                   start_codon_dict,
                   stop_codon_dict,
                   cds_dict,
                   exon_dict, db_name,
                   use_gene_name):
    """
    For each transcript, it extracts coordinates from all features (5'utr,
    3'utr, start codon, etc.) and merges and connects blocks of coordinates to
    create exon coordinates.
    """
    start_idx = 0
    stop_idx = 1
    
    component_list = [five_utr_dict, three_utr_dict, start_codon_dict,
                      stop_codon_dict, cds_dict, exon_dict]

    exon_table_rows = []
    for txt_id in transcript_id2strand:
        # blocks- list of coordinates of each component [(start, end)]
        # classification - the combined classification score for the transcript
        #                  as the highest number amongst all components.

        # For a given transcript, if it has CDS and exon features then exon
        # features take precedence
        if txt_id in exon_dict:
            this_component_list = component_list[0:-2]
            this_component_list.append(component_list[-1])
        elif txt_id in cds_dict:
            this_component_list = component_list[:-1]
        else:
            print "Transcript does not have a CDS or exon feature: %s" % txt_id
            sys.exit(1)

        exon_coords, classification = getAllComponents(txt_id,  
                                                  this_component_list)
        
#        exon_coords = mergeConnectBlocks(blocks) 
   
        # Getting gene Name from arbitrary CDS or exon entry 
        if txt_id in exon_dict:
            gene_name = getGeneName(exon_dict[txt_id][0], use_gene_name)
        else:
            gene_name = getGeneName(cds_dict[txt_id][0], use_gene_name)
            
        strand = transcript_id2strand[txt_id]

        for coord in exon_coords:
            if coord[start_idx] == 0:
                print "Error: GTF file has to be 1-based, not 0-based."
                sys.exit(1)
            exon_tuple = getExonTuple(db, txt_id, gene_name, chr, coord[start_idx],
                                    coord[stop_idx], strand, classification, db_name)
            exon_table_rows.append(exon_tuple)

    # Now insert all into the table
    db.insertListIntoTable("exon",
                           "transcript_id, gene_name, chr, start, end, strand, classification",
                           exon_table_rows,
                           db_name,
                           "?,?,?,?,?,?,?",
                           check_same_thread = False)

def createAnnotTables(db_obj, db_name):
    db_obj.createTable("exon",
        """exon_id INTEGER PRIMARY KEY AUTOINCREMENT, 
           transcript_id TEXT,
           gene_name TEXT,
           chr TEXT,
           start INTEGER,
           end INTEGER,
           strand TEXT,
           watson_strand_sequence TEXT,
           classification TEXT""",
         db_name,
         check_same_thread = False)

    db_obj.createIndex("exon_idx", "exon",
                       "transcript_id, gene_name", db_name,
                       check_same_thread = False)

    
    db_obj.createTable("intron",
        """intron_id INTEGER PRIMARY KEY AUTOINCREMENT,
           transcript_id TEXT,
           gene_name TEXT,
           chr TEXT,
           start INTEGER,
           end INTEGER,
           strand TEXT,
           watson_strand_sequence TEXT,
           classification TEXT""",
       db_name,
       check_same_thread = False)
    db_obj.createIndex("intron_idx", "intron",
                       "transcript_id, gene_name", db_name,
                       check_same_thread = False)

#   db_obj.createTable("cds",
#       """exon_id INTEGER PRIMARY KEY,
#          start INTEGER,
#          end INTEGER""",
#       db_name)

#   db_obj.createTable("start_codon",
#       """exon_id INTEGER PRIMARY KEY,
#          start INTEGER, 
#          end INTEGER""",
#       db_name)

#   db_obj.createTable("stop_codon",
#       """exon_id INTEGER PRIMARY KEY,
#          start INTEGER,
#          end INTEGER""",
#       db_name)

#   db_obj.createTable("five_utr",
#       """exon_id INTEGER PRIMARY KEY, 
#          start INTEGER,
#          end INTEGER""",
#       db_name)

#   db_obj.createTable("three_utr",
#       """exon_id INTEGER PRIMARY KEY,
#          start INTEGER,
#          end INTEGER""",
#       db_name)

    db_obj.createTable("gene",
        """gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
           name TEXT,
           chr TEXT,
           start INTEGER,
           end INTEGER,
           strand TEXT""",
        db_name,
        check_same_thread = False)
    db_obj.createIndex("gene_idx", "gene",
                       "name", db_name,
                        check_same_thread = False)

def buildGeneTable(db, db_name, chr):
    """
    This function derives the gene coordinates by finding the starting posiiton
    of the first exon(first occuring exon from all transcripts) and the ending
    position of the last exon(last occuring exon from all transcripts) for each
    gene.
    """
    # Local constants
    first_exon_idx = 0
    last_exon_idx = -1

    cg_nums = db.getDBRecords_Dict("""SELECT DISTINCT gene_name FROM exon 
                                      WHERE chr = \'%s\'""" % chr,
                                   db_name,
                                   check_same_thread = False)

    gene_rows = []
    for row in cg_nums:
        gene = row["gene_name"]


        # Need to escape names with "'"        
        mysql_query_gene = gene
        if "'" in mysql_query_gene:
            mysql_query_gene = mysql_query_gene.replace("'","''")

        # Get all the exons in this gene
        start_select_state = """SELECT start, strand, chr
                                FROM exon
                                WHERE gene_name=\'%s\'
                                ORDER BY start""" % mysql_query_gene

        start_records = db.getDBRecords_Dict(start_select_state, db_name,
                                             check_same_thread = False)

        # Get the strand of the gene, inferred from exon strands
        strand = getStrand(start_records, gene)

        if strand == None:
            continue

        end_select_state = """SELECT end FROM exon
                              WHERE gene_name=\'%s\'
                              ORDER BY end""" % mysql_query_gene
        end_records = db.getDBRecords_Dict(end_select_state, db_name,
                                           check_same_thread = False)

        gene_start = int(start_records[first_exon_idx]["start"])
        gene_end = int(end_records[last_exon_idx]["end"])

        # Insert gene information
        vals = (mysql_query_gene, strand, start_records[first_exon_idx]["chr"],
                gene_start, gene_end)
        gene_rows.append(vals)
        

    # Now insert all into the table
    db.insertListIntoTable("gene",
                           "name, strand, chr, start, end",
                           gene_rows,
                           db_name,
                           "?,?,?,?,?",
                           check_same_thread = False)

def extractInsertSeqs(tbl_name, db_obj, genome_file_name, db_name):
    """
    For each chromosome, gets all exon or intron coordinates and extracts the
    sequence.  The sequence is only extracted from the forward strand (watson
    strand).  No reverse complement operations are performed.
    """
    seq_file = open(genome_file_name)

    for tempSeqObj in SeqIO.parse(seq_file, "fasta"):

        # Get the chr name and sequence
        seqTitle = tempSeqObj.id.split()
        name = seqTitle[0]
        # Remove >
        chr_name = name.lstrip(">")

        scaff_seq = DNA(tempSeqObj.seq.tostring())

        # Get all exons or introns from this chr
        select_statement = """SELECT * FROM %s
                              WHERE chr=\'%s\'""" % (tbl_name,
                                                     chr_name)

        records = db_obj.getDBRecords_Dict(select_statement, db_name,
                                           check_same_thread = False)

        if records is None:
            print "No records found on %s" % chr_name
            continue

        for row in records:
            start = int(row["start"])
            end = int(row["end"])

            seq = scaff_seq.extractSeq(start, end)

            if tbl_name == "exon":
                update_statement = """UPDATE exon
                                      SET watson_strand_sequence=\'%s\'
                                      WHERE exon_id=%d""" % (seq,
                                                            int(row["exon_id"]))
            else: # tbl_name == intron
                update_statement = """UPDATE intron
                                      SET watson_strand_sequence=\'%s\'
                                      WHERE intron_id=%d""" % (seq,
                                                            int(row["intron_id"]))

            db_obj.updateTable(tbl_name, update_statement, db_name,
                               check_same_thread = False)
                

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getAllComponents(txt_id, dict_list):
    """
    Returns a list of coordinates (start, stop, classification) with the gff
    score, which in the case of the DMG2 annotations, is the classification
    I'm calling the coordinates blocks.
    """
    blocks = []

    classification = "."
    
    for i in range(len(dict_list)):
        if txt_id in dict_list[i]:
            for gff_obj in dict_list[i][txt_id]:
                blocks.append((int(gff_obj.start), 
                               int(gff_obj.end)))
                if type(gff_obj.score) == type(1):
                    if gff_obj.score != classification:
                        if classification == ".":
                            classification = gff_obj.score
                        elif int(gff_obj.score) < int(classification):
                            classification = gff_obj.score

    return blocks, classification

def getGeneName(gff_obj, use_gene_name):
    attributes = gff_obj.attributes.split(";")
    name = None
    if use_gene_name:
        for attrib in attributes:
            if "gene_name" in attrib:
                name = attrib
    if not name:
        name = attributes[0]
    name = name.split()[1]
    name = name.replace("\"","")
    return name

def getGTFLines(gtf_file):
    """
    Returns a dictionary associating chromosome to the gtf line
    {chr:[lines,]}
    """
    file = open(gtf_file)

    chr2gtf_lines = {}

    for line in file:
        # Remove comments
        if line.startswith("#"):
            continue
    
        line = formatLine(line)

        chr = line.split("\t")[0] 

        updateDictOfLists(chr2gtf_lines, chr, line)

    file.close()

    return chr2gtf_lines

def getStrand(db_records, gene_name):
    strandIsMinus = False
    strandIsPlus = False

    for row in db_records:
        if row["strand"] == "":
            continue
        if row["strand"] == "-":
            strandIsMinus = True
        if row["strand"] == "+":
            strandIsPlus = True

    # Error check
    if strandIsMinus and strandIsPlus:
        print "Disagreeing strands in gene: %s" % gene_name
        return None

    if strandIsMinus:
        return "-"
    if strandIsPlus:
        return "+"

    # Default return value
    return None


def get_transcript_id(gff_obj):
    t_id = gff_obj.attributes.split(";")[1]
    t_id = t_id.split()[1]
    t_id = t_id.replace("\"","")
    return t_id


def inferIntrons(db, db_name, chr):
    """
    Iterates through all transcripts, collects all exons, sorts them in order,
    and spaces in between each exon are considered intron coordinates.
    """

    # Get all transcripts
    txt_select = """SELECT DISTINCT transcript_id FROM exon
                    WHERE chr=\'%s\'""" % chr

    txt_records = db.getDBRecords_Dict(txt_select, db_name, check_same_thread = False)

    intron_rows = []
    for row in txt_records:
        txt = row["transcript_id"]

        # Escaping "'" in names
        if "'" in txt:
            txt = txt.replace("'","''") 

        # Get all exons of this transcript
        exon_select = """SELECT chr, start, end, strand, gene_name,
                                classification
                         FROM exon
                         WHERE transcript_id=\'%s\'
                         ORDER BY start""" % txt

        exon_records = db.getDBRecords_Dict(exon_select, db_name,
                                            check_same_thread = False)

        if not validTranscript(exon_records):
            print "Transcript has overlapping exons: %s" % txt
            continue

        prev_exon_end = 0
        for exon_row in exon_records:
            if prev_exon_end == 0:
                prev_exon_end = int(exon_row["end"])
                continue
            this_exon_start = int(exon_row["start"])

            intron_start = prev_exon_end + 1
            intron_end = this_exon_start - 1

            gene_name = exon_row["gene_name"]
            if "'" in gene_name:
                gene_name = gene_name.replace("'","''")

            vals = (txt, intron_start, intron_end, exon_row["chr"],
                    exon_row["strand"], gene_name, exon_row["classification"])

            intron_rows.append(vals)
    
            prev_exon_end = int(exon_row["end"]) 

    # Now insert all into the table
    db.insertListIntoTable("intron",
                           "transcript_id, start, end, chr, strand, gene_name, classification",
                           intron_rows,
                           db_name,
                           "?,?,?,?,?,?,?",
                           check_same_thread = False)

def getExonTuple(db, txt_id, gene_name, chr, start, end, strand, classification,
                db_name):

#   cols = """transcript_id, gene_name, chr, start, end, strand,
#             classification"""

    # If an apostrophe is in a name, then it needs to be escaped.
    if "'" in txt_id:
        txt_id = txt_id.replace("'","''")
    if "'" in gene_name:
        gene_name = gene_name.replace("'","''")

    vals = (txt_id, gene_name, chr, start, end, strand, classification)

    return vals

def insertFeature(db, table, table_dict, db_name):
    """
    Will insert the cds, 5UTR, 3UTR, start_codon, and stop_codon
    into the database.
    """ 
    for txt_id in table_dict:
        for gff_obj in table_dict[txt_id]:

            this_txt_id = get_transcript_id(gff_obj)
            if "'" in this_txt_id:
                this_txt_id = this_txt_id.replace("'","''")
            
            select_statement = """SELECT exon_id FROM exon
                                  WHERE
                                  transcript_id=\'%s\'
                                  AND
                                  start<=%d
                                  AND 
                                  end>=%d""" % (this_txt_id,
                                                int(gff_obj.start),
                                                int(gff_obj.end))
            record = db.getDBRow_Dict(select_statement, db_name,
                                      check_same_thread = False)
            if record == None:
                print "Could not find exon for %s:" % table
                print gff_obj
                continue

            exon_id = int(record["exon_id"])

            # Error checking before insertion
            row_exist_statement = """SELECT exon_id FROM %s
                                     WHERE exon_id=%d""" % (table, exon_id)
            if rowExists(db, db_name, row_exist_statement):
                print "This exon_id has multiple associated %s: %d" % (table,
                                                                       exon_id),
                print " in %s" % db_name
                continue

            cols = "exon_id, start, end"
            vals = "%d,%d,%d" % (exon_id, int(gff_obj.start), int(gff_obj.end))
            db.insertIntoTable(table, cols, vals, db_name, check_same_thread = False)

def mergeConnectBlocks(blocks):
    """
    Merges overlapping blocks and connects adjacent blocks.  This is done by
    creating an array that flags every position in the range of the blocks as
    containing a block portion.  Contiguous block flags are considered exonic
    regions.
    """
    start_idx = 0
    end_idx = 1

    # This dictionary is in the form {coord:0 (arbitrary value)}
    # The dictionary format prevents duplicates.
    exonic_dict = {} 
    
    exon_blocks = [] # Return value

    for block in blocks:
        for i in range(block[start_idx], block[end_idx] + 1):
            exonic_dict[i] = 0

    # Output exon locations
    exonic_locations = exonic_dict.keys() 
    exonic_locations.sort()

    lastPos = -1
    startPos = None

    for i in range(len(exonic_locations)):
        thisPos = exonic_locations[i]
        
        if thisPos != (lastPos + 1):
            if lastPos == -1:
                startPos = thisPos
                lastPos = thisPos
                continue
            # End of a block
            exon_blocks.append((startPos, lastPos))

            # update startPos for new block
            startPos = thisPos

        lastPos = thisPos

    # Insert last block
    exon_blocks.append((startPos, lastPos))

    return exon_blocks
        
def rowExists( db, db_name, select_statement):

    row = db.getDBRow_Dict( select_statement, db_name ,
                            check_same_thread = False)

    if row == None:
        return False
    else:
        return True

def validTranscript(exon_records):
    """
    Exon records should already be in order
    """
    last_start = -1
    last_end = -1

    for exon_row in exon_records:
        this_start = int(exon_row["start"])
        this_end = int(exon_row["end"])

        if this_start > this_end:
            return False

        if this_start > last_start and this_start > last_end:
            last_start = this_start
            last_end = this_end
        else:
            return False

    return True

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
