#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: heapSorter.py,v $
#   $Revision: 1.11 $
#
#   $Id: heapSorter.py,v 1.11 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the HeapSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter

#{
class HeapSorter(Sorter):
    """
    Heap sorter.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (HeapSorter) -> None
        Constructor.
        """
        super(HeapSorter, self).__init__()


    def percolateDown(self, i, length):
        """
        (HeapSorter, int, int) -> None
        Percolates the elements in the array with the given length
        and starting at the given position.
        """
        while 2 * i <= length:
            j = 2 * i
            if j < length and self._array[j + 1] \
                    > self._array[j]:
                j = j + 1
            if self._array[i] \
                    >= self._array[j]:
                break
            self.swap(i, j)
            i = j
#}>a

#{
    def buildHeap(self):
        """
        (HeapSorter) -> None
        Builds the heap.
        """
        i = self._n / 2
        while i > 0:
            self.percolateDown(i, self._n)
            i -= 1
#}>b

#{
    def _sort(self):
        """
        (HeapSorter) -> None
        Sorts the elements of the array.
        """
	base = self._array.baseIndex
	self._array.baseIndex = 1
        self.buildHeap()
        i = self._n
        while i >= 2:
            self.swap(i, 1)
            self.percolateDown(1, i - 1)
            i -= 1
	self._array.baseIndex = base
#}>c

    @staticmethod
    def main(*argv):
        "HeapSorter test program."
        print HeapSorter.main.__doc__
        Sorter.test(HeapSorter(), 100, 123)
        return 0

if __name__ == "__main__":
    sys.exit(HeapSorter.main(*sys.argv))
