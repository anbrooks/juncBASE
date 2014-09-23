#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: openScatterTableV2.py,v $
#   $Revision: 1.17 $
#
#   $Id: openScatterTableV2.py,v 1.17 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the OpenScatterTableV2 class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.17 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.hashTable import HashTable
from opus7.openScatterTable import OpenScatterTable
from opus7.exception import *

#{
class OpenScatterTableV2(OpenScatterTable):
    """
    Hash Table implemented as an open scatter table using an array.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, length = 0):
        """
        (OpenScatterTableV2 [, int]) -> None
        Constructs an open scatter table with the given length.
        """
        super(OpenScatterTableV2, self).__init__(length)

#{
    def withdraw(self, obj):
        """
        (OpenScatterTableV2, Object) -> None
        Withdraws the given object from this open scatter table.
        """
        if self._count == 0:
            raise ContainerEmpty
        i = self.findInstance(obj)
        if i < 0:
            raise KeyError
        while True:
            j = (i + 1) % len(self)
            while self._array[j]._state == self.OCCUPIED:
                h = self.h(self._array[j]._obj)
                if ((h <= i and i < j) or (i < j and j < h) or \
                        (j < h and h <= i)):
                    break
                j = (j + 1) % len(self)
            if self._array[j]._state == self.EMPTY:
                break
            self._array[i] = self._array[j]
            i = j
        self._array[i] = self.Entry(self.EMPTY, None)
        self._count -= 1
#}>a

    @staticmethod
    def main(*argv):
        "OpenScatterTableV2 test program."
        print OpenScatterTableV2.main.__doc__
        hashTable = OpenScatterTableV2(57)
        HashTable.test(hashTable)
        return 0

if __name__ == "__main__":
    sys.exit(OpenScatterTableV2.main(*sys.argv))
