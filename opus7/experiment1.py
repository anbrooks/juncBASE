#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: experiment1.py,v $
#   $Revision: 1.11 $
#
#   $Id: experiment1.py,v 1.11 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Experiment1 class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import example
from opus7.timer import Timer

class Experiment1(object):
    """
    Program that measures the running times of both
    a recursive and a non-recursive method to compute
    the Fibonacci numbers.
    """

    @staticmethod
    def main(*argv):
        "Experiment1 test program."
        print "3"
        print "n"
        print "fib1 s"
        print "fib2 s"
        timer1 = Timer()
        timer2 = Timer()
        for i in xrange(49):
            timer1.start()
            result = example.Fibonacci1(i)
            timer1.stop()

            timer2.start()
            result = example.Fibonacci2(i)
            timer2.stop()

            datum = "%d\t%g\t%g" % (
                i, timer1.getElapsedTime(), timer2.getElapsedTime())
            print datum
            sys.stderr.write(datum + "\n")
        return 0

if __name__== "__main__":
    sys.exit(Experiment1.main(*sys.argv))
