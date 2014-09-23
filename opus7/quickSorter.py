#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: quickSorter.py,v $
#   $Revision: 1.13 $
#
#   $Id: quickSorter.py,v 1.13 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the QuickSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.13 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.sorter import Sorter
from opus7.straightInsertionSorter import StraightInsertionSorter

#{
class QuickSorter(Sorter):
    """
    Base class from which all QuickSorter classes are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (QuickSorter) -> None
        Constructor.
        """
        super(QuickSorter, self).__init__()

    @abstractmethod
    def selectPivot(self): pass
#}>a

#{
    CUTOFF = 2 # minimum cut-off

    def quicksort(self, left, right):
        """
        (QuickSorter, left, right) -> None
        Recursively sorts the elements of the array between left and right.
        """
        if right - left + 1 > self.CUTOFF:
            p = self.selectPivot(left, right)
            self.swap(p, right)
            pivot = self._array[right]
            i = left
            j = right - 1
            while True:
                while i < j and self._array[i] < pivot:
                    i += 1
                while i < j and self._array[j] > pivot:
                    j -= 1
                if i >= j:
                    break
                self.swap(i, j)
                i += 1
                j -= 1
            if self._array[i] > pivot:
                self.swap(i, right)
            if left < i:
                self.quicksort(left, i - 1)
            if right > i:
                self.quicksort(i + 1, right)
#}>b

#{
    def _sort(self):
        """
        (QuickSorter) -> None
        Sorts the elements of the array.
        """
        self.quicksort(0, self._n - 1)
        sorter = StraightInsertionSorter()
        sorter.sort(self._array)
#}>c
