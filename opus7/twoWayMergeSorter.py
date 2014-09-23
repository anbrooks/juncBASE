#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: twoWayMergeSorter.py,v $
#   $Revision: 1.10 $
#
#   $Id: twoWayMergeSorter.py,v 1.10 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the TwoWayMergeSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter
from opus7.array import Array

#{
class TwoWayMergeSorter(Sorter):
    """
    Two-way merge sorter.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (TwoWayMergeSorter) -> None
        Constructor.
        """
        super(TwoWayMergeSorter, self).__init__()
        self._tempArray = None
#}>a

#{
    def merge(self, left, middle, right):
        """
        (TwoWayMergeSorter, int, int, int) -> None
        Merges two sorted sub-arrays,
        array[left] ... array[middle] and
        array[middle + 1] ... array[right]
        using the temporary array.
        """
        i = left
        j = left
        k = middle + 1
        while j <= middle and k <= right:
            if self._array[j] < self._array[k]:
                self._tempArray[i] = self._array[j]
                i += 1
                j += 1
            else:
                self._tempArray[i] = self._array[k]
                i += 1
                k += 1
        while j <= middle:
            self._tempArray[i] = self._array[j]
            i += 1
            j += 1
        for i in xrange(left, k):
            self._array[i] = self._tempArray[i]
#}>b

#{
    def _sort(self):
        """
        (TwoWayMergeSorter) -> None
        Sorts the elements of the array.
        """
        self._tempArray = Array(self._n)
        self.mergesort(0, self._n - 1)
        self._tempArray = None

    def mergesort(self, left, right):
        """
        (TwoWayMergeSorter, int, int) -> None
        Sorts the elements of the array array[left] ... array[right].
        """
        if left < right:
            middle = (left + right) / 2
            self.mergesort(left, middle)
            self.mergesort(middle + 1, right)
            self.merge(left, middle, right)
#}>c

    @staticmethod
    def main(*argv):
        "TwoWayMergeSorter test program."
        print TwoWayMergeSorter.main.__doc__
        Sorter.test(TwoWayMergeSorter(), 100, 123)
        return 0

if __name__ == "__main__":
    sys.exit(TwoWayMergeSorter.main(*sys.argv))
