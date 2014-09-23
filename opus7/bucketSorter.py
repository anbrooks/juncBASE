#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: bucketSorter.py,v $
#   $Revision: 1.10 $
#
#   $Id: bucketSorter.py,v 1.10 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BucketSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter
from opus7.array import Array

#{
class BucketSorter(Sorter):
    """
    Bucket sorter.
    """

#}@head

#{

    # ...
#!}
#}@tail

#{
    def __init__(self, m):
        """
        (BucketSorter, int) -> None
        Constructs a bucket sorter with the given number of buckets.
        """
        super(BucketSorter, self).__init__()
        self._m = m
        self._count = Array(self._m)
#}>a

#{
    def _sort(self):
        """
        (BucketSorter) -> None
        Sorts the elements of the array.
        """
        for i in xrange(self._m):
            self._count[i] = 0
        for j in xrange(self._n):
            self._count[self._array[j]] += 1
        j = 0
        for i in xrange(self._m):
            while self._count[i] > 0:
                self._array[j] = i
                j += 1
                self._count[i] -= 1
#}>b
    
    @staticmethod
    def main(*argv):
        "BucketSorter test program."
        print BucketSorter.main.__doc__
        Sorter.test(BucketSorter(1024), 100, 123, 1024)
        return 0

if __name__ == "__main__":
    sys.exit(BucketSorter.main(*sys.argv))
