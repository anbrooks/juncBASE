#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: bubbleSorter.py,v $
#   $Revision: 1.10 $
#
#   $Id: bubbleSorter.py,v 1.10 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BubbleSorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sorter import Sorter

#{
class BubbleSorter(Sorter):
    """
    Bubble sorter.
    """

    def __init__(self):
        """
        (BubbleSorter) -> None
        Constructor.
        """
        super(BubbleSorter, self).__init__()

    def _sort(self):
        """
        (BubbleSorter) -> None
        Sorts the elements of the array.
        """
        i = self._n
        while i > 1:
            for j in xrange(i - 1):
                if self._array[j] > self._array[j + 1]:
                    self.swap(j, j + 1)
            i -= 1
#}>a

    @staticmethod
    def main(*argv):
        "BubbleSorter test program."
        print BubbleSorter.main.__doc__
        Sorter.test(BubbleSorter(), 100, 123)
        return 0

if __name__ == "__main__":
    sys.exit(BubbleSorter.main(*sys.argv))
