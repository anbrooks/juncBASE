#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: experiment2.py,v $
#   $Revision: 1.11 $
#
#   $Id: experiment2.py,v 1.11 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Experiment2 class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.demo9 import Demo9

class Experiment2(object):
    """
    Program that measures the running times of various sorting algorithms.
    """

    @staticmethod
    def main(*argv):
        "Experiment2 test program."
        print "4"
        print "sort"
        print "length"
        print "seed"
        print "time"
        seeds = ["1", "57", "12345", "7252795", "3127"]
        for seed in seeds:
            lengths = ["10", "25", "50", "75",
                "100", "250", "500", "750",
                "1000", "1250", "1500", "1750", "2000"]
            for length in lengths:
                Demo9.main("Demo9.main", length, seed, "7")
            lengths = ["3000", "4000", "5000", "6000",
                "7000", "8000", "9000", "10000"]
            for length in lengths:
                Demo9.main("Demo9.main", length, seed, "3")
            lengths = ["20000", "30000", "40000", "50000",
                "60000", "70000", "80000", "90000", "100000"]
            for length in lengths:
                Demo9.main("Demo9.main", length, seed, "1")
        return 0

if __name__== "__main__":
    sys.exit(Experiment2.main(*sys.argv))
