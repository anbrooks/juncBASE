# pysqlite_wrap.py
# Author: Angela Brooks
# Program Completion Date:
# Description: This is a wrapper class for the pysqlite module.  I noticed that I
# us the same sort of functions with the module, so this will be a faster way
# to interface with a SQLite database
# Modification Date(s):
#
# class DB:
# 	__init__(self,host,user,passwd)
#	createDatabase(self,database_name)
#	createTable(self, table_name, create_definition, database)
#	insertIntoTable(self, table_name, column_names, values, database )
#	getDBRecords_Dict( self, select_statement, database)
#	getDBRow_Dict( self, select_statement, database)
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

import sys, getopt, pdb, os

from pysqlite2 import dbapi2 as sqlite
###########
#CONSTANTS#
###########

    
###############
#END CONSTANTS#
###############


#########
#CLASSES#
#########
def dict_factory(cursor, row):
		d = {}
		for idx, col in enumerate(cursor.description):
			d[col[0]] = row[idx]
		return d
class DictCursor(sqlite.Cursor):
	def __init__(self, *args, **kwargs):
		sqlite.Cursor.__init__(self, *args, **kwargs)
		self.row_factory = dict_factory

class DB:

	def __init__(self, db_dir="./"):
		self.db_dir=db_dir
		if not db_dir.endswith("/"):
			self.db_dir=db_dir + "/"

	def __connect(self, database, check_same_thread = True):
		return sqlite.connect( self.db_dir + database, check_same_thread )	

	def createDatabase(self,database_name):
		# Delete previoius database
		if os.path.exists(self.db_dir + database_name):
			os.remove(self.db_dir + database_name)
			
		conn = sqlite.connect( self.db_dir + database_name )
		cursor = conn.cursor()
	
		conn.commit()
		cursor.close()
		conn.close()
		
	def createIndex(self,index_name, table_name, comma_separated_cols,
					database_name, check_same_thread = True):
			
		conn = self.__connect( database_name, check_same_thread )
		cursor = conn.cursor()
	
		create_index = "CREATE INDEX %s ON %s " % (index_name, table_name)
		create_index += "(%s)" % comma_separated_cols
		
		cursor.execute( create_index )
	
		conn.commit()
		cursor.close()
		conn.close()

	def createTable(self, table_name, create_definition, database,
                    check_same_thread = True):
		conn = self.__connect( database, check_same_thread )
		cursor = conn.cursor()
	
		#Delete previous table
		try:
			drop_statement = "DROP TABLE IF EXISTS %s" % table_name
			cursor.execute( drop_statement )
			conn.commit()

		except sqlite.IntegrityError, e:
			print "Create Table Error %d: %s" % (e.args[0], e.args[1])

		create_statement = "CREATE TABLE %s (%s)" % (table_name,
							create_definition)
		cursor.execute( create_statement )
		conn.commit()
		cursor.close()
		conn.close()
	
	def insertIntoTable(self, table_name, column_names, values, database,
                      val_string=None, check_same_thread = True ):
		conn = self.__connect( database, check_same_thread )
		cursor = conn.cursor()
		
		if not val_string:
			val_string_comp = ["?" for n in range(values.count(",")+1)]
			val_string = ",".join(val_string_comp)

		insert_statement = "INSERT INTO %s (%s) VALUES (%s)" % (table_name,
                                                            column_names,
								  	                                        val_string)
		try:
			cursor.execute( insert_statement, values )
			conn.commit()
		except sqlite.IntegrityError, e:
			print "Insert Error %d: %s" % (e.args[0], e.args[1])
			
		cursor.close()
		conn.close()	

	def insertListIntoTable(self, table_name, column_names, values_list, database,
                          val_string=None, check_same_thread = True):
		conn = self.__connect(database, check_same_thread)
		cursor = conn.cursor()

		if not val_string:
			val_string_comp = ["?" for n in range(values.count(",")+1)]
			val_string = ",".join(val_string_comp)

		insert_statement = "INSERT INTO %s (%s) VALUES (%s)" % (table_name,
                                                            column_names,
														                                val_string)
		try:
			cursor.executemany( insert_statement, values_list )
			conn.commit()
		except sqlite.IntegrityError, e:
			print "Insert Error %d: %s" % (e.args[0], e.args[1])
			
		cursor.close()
		conn.close()	

	def updateTable(self, table_name, update_statement, database,
                    check_same_thread = True):
		conn = self.__connect(database, check_same_thread)
		cursor = conn.cursor()

		try:
			cursor.execute(update_statement)
			conn.commit()
		except sqlite.IntegrityError, e:
			print "Update Error %d: %s" % (e.args[0], e.args[1])

		cursor.close()
		conn.close()
	
	def getDBRecords_Dict( self, select_statement, database, check_same_thread = True):
		"""
		Returns [] if no records exist.
		"""
        	#Connect to Database
        	conn = self.__connect(database, check_same_thread)

        	conn.row_factory = sqlite.Row
        	cursor = conn.cursor()
#	        cursor = conn.cursor(factory=DictCursor)
		try:
			cursor.execute ( select_statement )
			return_result =  cursor.fetchall()
		
		except sqlite.IntegrityError, e:
			print "Get Records Error %d: %s" % (e.args[0], e.args[1])
			return None

        	#Close connections
        	cursor.close()
        	conn.close()

        	return return_result

	def getDBRow_Dict( self, select_statement, database, check_same_thread = True):
		"""
		Returns None if no row exists
		"""
                #Connect to Database
		conn = self.__connect(database, check_same_thread)

		conn.row_factory = sqlite.Row
		cursor = conn.cursor()
#		cursor = conn.cursor(factory=DictCursor)

		try:
			cursor.execute ( select_statement )
			return_result =  cursor.fetchone()
		except sqlite.IntegrityError, e:
			print "Get Row Error %d: %s" % (e.args[0], e.args[1])
			return None
		
                cursor.close()
                conn.close()
                                                                                
                return return_result
	

#############
#END CLASSES#
#############
 
######
#MAIN#	
######
def main():		

	db = DB()

	db.createDatabase("lotsOstuff")

#	db.createTable( "pens", "id INTEGER PRIMARY KEY AUTOINCREMENT, color VARCHAR(10)", "lotsOstuff")
	db.createTable( "pens", "color VARCHAR(10)", "lotsOstuff")
	db.insertIntoTable( "pens", "color",('blue',), "lotsOstuff", "?")
	db.insertIntoTable( "pens", "color", ('red',), "lotsOstuff")

	db.insertListIntoTable("pens", "color", [('pink',), 
                                           ('green',),
                                           ('black',),
                                           ('yellow',),
                                           ('magenta',),
                                           ('silver',),
                                           ('gold',),
                                           ('purple',)],
                          "lotsOstuff",
                          "?")
	select_statement = "SELECT * FROM pens WHERE color=\'blue\'"	

	pen_record = db.getDBRow_Dict( select_statement, "lotsOstuff" )

#	print "id=%d" % int(pen_record["id"])
	print "color=%s" % pen_record["color"]

	select_statement2 = "SELECT * FROM pens"
	
	pen_record = db.getDBRecords_Dict( select_statement2, "lotsOstuff")

	print "The pens are:"
	print pen_record

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
