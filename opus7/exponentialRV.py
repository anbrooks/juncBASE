#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: exponentialRV.py,v $
#   $Revision: 1.9 $
#
#   $Id: exponentialRV.py,v 1.9 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the ExponentialRV class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.9 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import math
from opus7.randomVariable import RandomVariable
from opus7.randomNumberGenerator import RandomNumberGenerator

#{
class ExponentialRV(RandomVariable):
    """
    Exponentially distributed random variable.
    """

    def __init__(self, mu):
        """
        (ExponentialRV, double) -> None
        Constructs an exponentially distributed random variable
        with the given mean.
        """
        self._mu = mu

    def getNext(self):
        """
        (ExponentialRV) -> double
        Returns the next sample.
        """
        return -self._mu * math.log(RandomNumberGenerator.next)
#}>a

    @staticmethod
    def main(*argv):
	"ExponentialRV test program."
	print ExponentialRV.main.__doc__
	rv = ExponentialRV(100)
	for i in range(10):
	    print rv.next
	return 0

if __name__== "__main__":
    sys.exit(ExponentialRV.main(*sys.argv))
