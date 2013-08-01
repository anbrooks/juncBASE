#!/lab/64/bin/python
# helperFunctions.py
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import random
import math
import os
import time

from subprocess import Popen, PIPE
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

class CommentedFile:
    """
    Taken from:
    http://www.mfasold.net/blog/2010/02/python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines/
    """
    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring):
            line = self.f.next()
        return line
    def __iter__(self):
        return self
    def close(self):
        self.f.close()

###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
#    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
#    opt_parser.add_option("<opt>",
#                          dest=,
#                          type=,
#                          help=,)

#    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
#    opt_parser.check_required("<req_arg>")
			
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

def coordsOverlap(start1, end1, start2, end2):
    """
    Just check if the coordinates overlap.
    """
    if start1 >= start2 and start1 <= end2:
        return True
    if end1 >= start2 and end1 <= end2:
        return True
    if start1 <= start2 and end1 >= end2:
        return True

    return False

def median(vals):
    vals.sort()
    val_len = len(vals)

    if val_len % 2 != 0:
        return vals[(val_len + 1)/2 - 1]
    else:
        lower = vals[val_len/2-1]
        upper = vals[val_len/2]
        return float(lower + upper)/2
        
def runCmd(cmd_str, which_shell, pCommand=False):
    try:
        if pCommand:
            print cmd_str

        p = Popen(cmd_str, shell=True, executable=which_shell)
        sts = os.waitpid(p.pid,0)
    except:
        print "Something wrong with: %s" % cmd_str
        sys.exit(1)
        
def getShannonIndex(pos2countDict, totalCount):
    """
    Calculates a Shannon-Wiener Index for the positions and counts.  It treats
    each position like a unique species.
    """

    if totalCount == 0:
        return 0

    summation = 0

    for pos in pos2countDict:
        p = float(pos2countDict[pos]) / totalCount

        if p == 0.0:
            continue

        summation += p * math.log(p, 2)

    if summation == 0:
        return 0

    return -summation

def runLSF(cmd, bsub_out, job_id, queue="week", mem=8, group="cancerfolk",
tmp_file_name = "tmp.bjob", job_group=None):
    bsub_options = "#!/bin/tcsh\n"
    bsub_options += "#BSUB -q %s\n" % queue
    if queue == "hour":
        bsub_options += "#BSUB -W 240\n"
    bsub_options += "#BSUB -R \"rusage[mem=%d]\"\n" % mem
    bsub_options += "#BSUB -P %s\n" % group

    bsub_options += "#BSUB -o %s\n" % bsub_out
    bsub_options += "#BSUB -J %s\n" % job_id

    if job_group:
        bsub_options += "#BSUB -g /%s\n" % job_group

    # Make tmp file
    tmp_file = open(tmp_file_name, "w")
        
    tmp_file.write(bsub_options)

    tmp_file.write(cmd)
     
    tmp_file.close()

    os.system("bsub < %s" % tmp_file_name)

    # Remove tmp_file
    os.remove(tmp_file_name)

def updateDictOfLists(d, key, item):
    """
    I noticed that in many scripts I have a dictionary of lists.  When I want
    to add to the dictionary, I always have to first check if the key exists,
    if it doesn't I have to map the new key to a list, but if it does, I just
    append to the existing list.  Doing this many times in separate pieces of
    code makes this function useful.
    """ 
    try:
        d[key].append(item)
    except KeyError:
        d[key] = [item]

def updateDictOfLists_extend(d, key, item):
    """
    I noticed that in many scripts I have a dictionary of lists.  When I want
    to add to the dictionary, I always have to first check if the key exists,
    if it doesn't I have to map the new key to a list, but if it does, I just
    append to the existing list.  Doing this many times in separate pieces of
    code makes this function useful.

    In some cases, the item might be a list itself, which I want to extend to
    the existing list instead of just append a list to a list.
    This function would also work if you wanted to just add one item.
    """ 
    try:
        if type(item) == type([]):
            d[key].extend(item)
        else:
            d[key].append(item)
    except KeyError:
        if type(item) == type([]):
            d[key] = item 
        else:
            d[key] = [item]

def updateDictOfSets(d, key, item):
    """
    Similar functionality as updateDictOfLists
    """
    try:
        d[key].add(item)
    except KeyError:
        d[key] = set([item])

def waitForChildren(children_processes, sleeptime=5):
    ctr = 0
    num_children = len(children_processes)
    while 1:
        for child in children_processes:
            if not child.poll() is None:
                ctr += 1
            sys.stdout.flush()

        if ctr == num_children:
            break

        time.sleep(sleeptime)
        

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
