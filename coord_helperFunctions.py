#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
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

from helperFunctions import coordsOverlap

from opus7.avlTree import AVLTree
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
    opt_parser.add_option("-r",
                          dest="run_test",
                          action="store_true",
                          help="Run tests of functions",
                          default=None)

    (options, args) = opt_parser.parse_args()

    coord_list = [(10,20),
                  (22, 25),
                  (80, 100),
                  (200, 280)]

    coord_tree = getSearchTree(coord_list)

    print "Test 1"
    if hasOuterContainer(coord_tree, (8,25)):
        print "Pass"
    else:
        print "Fail"

    print ""

    print "Test 2"
    if hasOuterContainer(coord_tree, (150,250)):
        print "Fail"
    else:
        print "Pass"
	
    print ""


    print "Test 3"
    if hasOuterContainer(coord_tree, (74,180)):
        print "Pass"
    else:
        print "Fail"

    print ""

    print "Test 4"
    if hasOuterContainer(coord_tree, (30,60)):
        print "Fail"
    else:
        print "Pass"

    print
    print "List of coords: %s" % coord_list

    print "Test 5"
    internalCoords1 = []   
    findInternalCoords(coord_tree, internalCoords1, (8,25))
    print "Coords in between (8,25):"
    print internalCoords1
    
    print "Test 6"
    internalCoords2 = []   
    findInternalCoords(coord_tree, internalCoords2, (50,75))
    print "Coords in between (50,75):"
    print internalCoords2

    print "Test 7"
    internalCoords3 = []   
    findInternalCoords(coord_tree, internalCoords3, (50,150))
    print "Coords in between (50,150):"
    print internalCoords3
    print

    print "Test 8"
    if hasOverlap(coord_tree, (8,25)):
        print "Pass"
    else:
        print "Fail"
    print

    print "Test 9"
    if hasOverlap(coord_tree, (50,75)):
        print "Fail"
    else:
        print "Pass"
    print

    print "Test 10"
    if hasOverlap(coord_tree, (70,200)):
        print "Pass"
    else:
        print "Fail"
    print

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def hasOuterContainer(coordTree, this_coord):

    this_tree = coordTree
    try:
        # Make sure that the tree is nonEmpty
        this_key = this_tree.getKey()
    except:
        return False

    if outerContained(this_coord, this_key):
        return True

    if not this_tree.getIsLeaf():
        if this_coord < this_key:
            if this_tree.getLeft().isEmpty:
                return False
            return hasOuterContainer(this_tree.getLeft(), this_coord)

        if this_coord > this_key:
            if this_tree.getRight().isEmpty:
                return False
            return hasOuterContainer(this_tree.getRight(), this_coord)

    return False

def hasOverlap(coordTree, this_coord):

    this_tree = coordTree
    try:
        # Make sure that the tree is nonEmpty
        tree_coord = this_tree.getKey()
    except:
        return False

    if coordsOverlap(this_coord[0], 
                     this_coord[1],
                     tree_coord[0],
                     tree_coord[1]):
        return True

    if not this_tree.getIsLeaf():
        if this_coord < this_tree.getKey():
            if this_tree.getLeft().isEmpty:
                return False
            return hasOverlap(this_tree.getLeft(), this_coord)

        if this_coord > this_tree.getKey():
            if this_tree.getRight().isEmpty:
                return False
            return hasOverlap(this_tree.getRight(), this_coord)

    return False

def findInternalCoords(coordTree, coord_list, this_coord):

    this_tree = coordTree
    try:
        this_key = this_tree.getKey()
    except:
        return

    thisIsContained = False
    if outerContained(this_coord, this_key):
        coord_list.append(this_tree.getKey())
        thisIsContained = True

    if not this_tree.getIsLeaf():
        if this_coord < this_key or thisIsContained:
            if this_tree.getLeft().isEmpty:
                return
            return findInternalCoords(this_tree.getLeft(), coord_list, this_coord)

        if this_coord > this_key or thisIsContained:
            if this_tree.getRight().isEmpty:
                return
            return findInternalCoords(this_tree.getRight(), coord_list, this_coord)

def findOverlapCoords(coordTree, coord_list, this_coord):

    this_tree = coordTree
    try:
        this_key = this_tree.getKey()
    except:
        return

    thisIsOverlap = False
    if coordsOverlap(this_coord[0],
                     this_coord[1],
                     this_key[0],
                     this_key[1]):
        coord_list.append(this_tree.getKey())
        thisIsOverlap = True

    if not this_tree.getIsLeaf():
        if this_coord < this_key or thisIsOverlap:
            if this_tree.getLeft().isEmpty:
                return
            return findOverlapCoords(this_tree.getLeft(), coord_list, this_coord)

        if this_coord > this_key or thisIsOverlap:
            if this_tree.getRight().isEmpty:
                return
            return findOverlapCoords(this_tree.getRight(), coord_list, this_coord)

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getSearchTree(coord_list):
    """
    coord_list is a list of (start, end) positions.
    """
    coordTree = AVLTree()

    # Make sure coords are unique
    coord_list_set = set(coord_list)
    
    for coord in coord_list_set:
        # ERROR CHECKING
        if type(coord) != type((1,2)):
            raise TypeError, ("List should be of tuples")
        if len(coord) != 2:
            raise ValueError, ("Tuple should only have two elements)")
        # Will raise exception if the coordinates are not numbers
        val_check = coord[1] - coord[0]

        coordTree.insert(coord)
       
    return coordTree 
    
def outerContained(this_coord, tree_coord):

    if this_coord[0] <= tree_coord[0]:
        if this_coord[1] >= tree_coord[1]:
            return True

    return False 

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
