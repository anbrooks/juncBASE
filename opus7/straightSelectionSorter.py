#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: straightSelectionSorter.py,v $
#   $Revision: 1.10 $
#
#   $Id: straightSelectionSorter.py,v 1.10 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the StraightSelectionSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter

#{
class StraightSelectionSorter(Sorter):
    """
    Straight selection sorter.
    """

    def __init__(self):
        """
        (StraightSelectionSorter) -> None
        Constructor.
        """
        super(StraightSelectionSorter, self).__init__()

    def _sort(self):
        """
        (StraightSelectionSorter) -> None
        Sorts the elements of the array.
        """
        i = self._n
        while i > 1:
            max = 0
            for j in xrange(i):
                if self._array[j] > self._array[max]:
                    max = j
            self.swap(i - 1, max)
            i -= 1
#}>a

    @staticmethod
    def main(*argv):
        "StraightSelectionSorter test program."
        print StraightSelectionSorter.main.__doc__
        Sorter.test(StraightSelectionSorter(), 100, 123)
        return 0

if __name__ == "__main__":
    sys.exit(StraightSelectionSorter.main(*sys.argv))
