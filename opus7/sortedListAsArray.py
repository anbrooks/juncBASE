#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sortedListAsArray.py,v $
#   $Revision: 1.19 $
#
#   $Id: sortedListAsArray.py,v 1.19 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SortedListAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.19 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sortedList import SortedList
from opus7.orderedListAsArray import OrderedListAsArray
from opus7.exception import *

#{
class SortedListAsArray(OrderedListAsArray, SortedList):
    """
    A sorted list implemented using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size = 0):
        """
        (SortedListAsArray [, int]) -> None
        Constructs a sorted list of the given size.
        """
        super(SortedListAsArray, self).__init__(size)
#}>a

#{
    def insert(self, obj):
        """
        (SortedListAsArray, Object) -> None
        Inserts the given object into this sorted list.
        """
        if self._count == len(self._array):
            raise ContainerFull
        i = self._count
        while i > 0 and self._array[i - 1] > obj:
            self._array[i] = self._array[i - 1]
            i -= 1
        self._array[i] = obj
        self._count += 1
#}>b

#{
    def findOffset(self, obj):
        """
        (SortedListAsArray, Object) -> int
        Finds the offset in this list
        of the object that equals the given object.
        """
        left = 0
        right = self._count - 1
        while left <= right:
            middle = (left + right) / 2
            if obj > self._array[middle]:
                left = middle + 1
            elif obj < self._array[middle]:
                right = middle - 1
            else:
                return middle
        return -1
#}>c

#{
    def find(self, obj):
        """
        (SortedListAsArray, Object) -> Object
        Finds the object in this list that equals the given object.
        """
        offset = self.findOffset(obj)
        if offset >= 0:
            return self._array[offset]
        else:
            return None

    class Cursor(OrderedListAsArray.Cursor):
        """
        A cursor that refers to an object in a sorted list.
        """

        def __init__(self, list, offset):
            """
            (SortedListAsArray.Cursor, SortedListAsArray, int) -> Object
            Constructs a cursor that refers to the object
            at the given offset in the given list.
            """
            super(SortedListAsArray.Cursor, self) \
                .__init__(list, offset)

        def insertAfter(self, obj):
            """
            (SortedListAsArray.Cursor, Object) -> None
            Not allowed in sorted lists.
            """
            raise TypeError

        def insertBefore(self, obj):
            """
            (SortedListAsArray.Cursor, Object) -> None
            Not allowed in sorted lists.
            """
            raise TypeError

    def findPosition(self, obj):
        """
        (SortedListAsArray, Object) -> SortedListAsArray.Cursor
        Finds the position of an object in this sorted list
        that equals the given object and returns a cursor
        that refers to that object.
        """
        return self.Cursor(self, self.findOffset(obj))
#}>d

#{
    def withdraw(self, obj):
        """
        (SortedListAsArray, Object) -> None
        Withdraws the given object from this sorted list.
        """
        if self._count == 0:
            raise ContainerEmpty
        offset = self.findOffset(obj)
        if offset < 0:
            raise KeyError
        i = offset
        while i < self._count:
            self._array[i] = self._array[i + 1]
            i += 1
        self._array[i] = None
        self._count -= 1
#}>e

    @staticmethod
    def main(*argv):
        "SortedListAsArray test program."
        print SortedListAsArray.main.__doc__
        slist = SortedListAsArray(10)
        SortedList.test(slist)
        return 0

if __name__ == "__main__":
    sys.exit(SortedListAsArray.main(*sys.argv))
