#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: radixSorter.py,v $
#   $Revision: 1.11 $
#
#   $Id: radixSorter.py,v 1.11 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the RadixSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter
from opus7.array import Array

#{
class RadixSorter(Sorter):
    """
    Radix sorter.
    """

#}@head

#{

    # ...
#}@tail

#{
    r = 8
    R = 1 << r
    p = (32 + r - 1) / r

    def __init__(self):
        """
        (RadixSorter) -> None
        Constructor.
        """
        self._count = Array(self.R)
        self._tempArray = None
#}>a

#{
    def _sort(self):
        """
        (RadixSorter) -> None
        Sorts the elements of the array.
        """
        self._tempArray = Array(self._n)
        for i in xrange(self.p):
            for j in xrange(self.R):
                self._count[j] = 0
            for k in xrange(self._n):
                self._count[(self._array[k] \
                    >> (self.r*i)) & (self.R-1)] += 1
                self._tempArray[k] = self._array[k]
            pos = 0
            for j in xrange(self.R):
                tmp = pos
                pos += self._count[j]
                self._count[j] = tmp
            for k in xrange(self._n):
                j = (self._tempArray[k] \
                    >> (self.r*i)) & (self.R-1)
                self._array[self._count[j]] = self._tempArray[k]
                self._count[j] += 1
#}>b

    @staticmethod
    def main(*argv):
        "RadixSorter test program."
        print RadixSorter.main.__doc__
        Sorter.test(RadixSorter(), 100, 123)
        return 0

if __name__ == "__main__":
    sys.exit(RadixSorter.main(*sys.argv))
