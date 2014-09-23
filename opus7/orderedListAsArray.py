#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: orderedListAsArray.py,v $
#   $Revision: 1.36 $
#
#   $Id: orderedListAsArray.py,v 1.36 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the OrderedListAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.36 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import exceptions
from opus7.orderedList import OrderedList
from opus7.array import Array
from opus7.cursor import Cursor
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class OrderedListAsArray(OrderedList):
    """
    Ordered list implemented using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size = 0):
        """
        (OrderedListAsArray [, int]) -> None
        Constructs an ordered list of the given size.
        """
        super(OrderedListAsArray, self).__init__()
        self._array = Array(size)
#}>a

#{
    def insert(self, obj):
        """
        (OrderedListAsArray, Object) -> None
        Inserts the given object at the end of this list.
        """
        if self._count == len(self._array):
            raise ContainerFull
        self._array[self._count] = obj
        self._count += 1
#}>b

    def purge(self):
        """
        (OrderedListAsArray) -> None
        Purges this ordered list.
        """
        while self._count > 0:
            self._count = self._count = 1
            self._array[self._count] = None

    def getIsFull(self):
        """
        (OrderedListAsArray) -> bool
        Returns true if this ordered list is full.
        """
        return self._count == len(self._array)

    def accept(self, visitor):
        """
        (OrderedListAsArray, Visitor) -> bool
        Makes the given visitor visit the objects in this ordered list.
        """
        assert isinstance(visitor, Visitor)
        for i in xrange(self._count):
            visitor.visit(self._array[i])
            if visitor.isDone:
                return
#{
    def __contains__(self, obj):
        """
        (OrderedListAsArray, Object) -> bool
        Returns true if the given object instance is in this ordered list.
        """
        for i in xrange(self._count):
            if self._array[i] is obj:
                return True
        return False

    def find(self, obj):
        """
        (OrderedListAsArray, Object) -> Object
        Finds an object in this ordered list that equals the given object.
        """
        for i in xrange(self._count):
            if self._array[i] == obj:
                return self._array[i]
        return None
#}>c

#{
    def withdraw(self, obj):
        """
        (OrderedListAsArray, Object) -> None
        Withdraws the given object instance from this ordered list.
        """
        if self._count == 0:
            raise ContainerEmpty
        i = 0
        while i < self._count and self._array[i] is not obj:
            i += 1
        if i == self._count:
            raise KeyError
        while i < self._count - 1:
            self._array[i] = self._array[i + 1]
            i += 1
        self._array[i] = None
        self._count -= 1
#}>d

#{
    def findPosition(self, obj):
        """
        (OrderedListAsArray, Object) -> OrderedListAsArray.Cursor
        Finds the position of an object in this list
        that equals the given object and returns a cursor
        that refers to that object.
        """
        i = 0
        while i < self._count and self._array[i] != obj:
            i += 1
        return self.Cursor(self, i)

    def __getitem__(self, offset):
        """
        (OrderedListAsArray, int) -> Object
        Returns the object in this list at the given._offset.
        """
        if offset < 0 or offset >= self._count:
            raise IndexError
        return self._array[offset]
#}>e

#{
    class Cursor(Cursor):
        """
        A cursor that refers to an object in an ordered list.
        """

#}@head

#{

        # ...
#}@tail

#{
        def __init__(self, list, offset):
            """
            (OrderedListAsArray.Cursor, OrderedListAsArray, int) -> None
            Constructs a cursor that refers to the object
            at the given._offset of the given list.
            """
	    super(OrderedListAsArray.Cursor, self).__init__(list)
            self._offset = offset

        def getDatum(self):
            """
            (OrderedListAsArray.Cursor) -> Object
            Returns the object to which this cursor refers.
            """
            if self._offset < 0 \
                    or self._offset >= self._list._count:
                raise IndexError
            return self._list._array[self._offset]
#}>h

#{
        def insertAfter(self, obj):
            """
            (OrderedListAsArray.Cursor, Object) -> None
            Inserts the given object into the list
            after the object to which this cursor refers.
            """
            if self._offset < 0 \
                    or self._offset >= self._list._count:
                raise IndexError
            if self._list._count == len(self._list._array):
                raise ContainerFull
            insertPosition = self._offset + 1
            i = self._list._count
            while i > insertPosition:
                self._list._array[i] = self._list._array[i - 1]
                i -= 1
            self._list._array[insertPosition] = obj
            self._list._count += 1
#}>f

        def insertBefore(self, obj):
            """
            (OrderedListAsArray.Cursor, Object) -> None
            Inserts the given object into the list
            before the object to which this cursor refers.
            """
            if self._offset < 0 \
                    or self._offset >= self._list._count:
                raise IndexError
            if self._list._count == len(self._list._array):
                raise ContainerFull
            insertPosition = self._offset
            i = self._list._count
            while i > insertPosition:
                self._list._array[i] = self._list._array[i - 1]
                i -= 1
            self._list._array[insertPosition] = obj
            self._list._count += 1
            self._offset += 1

#{
        def withdraw(self):
            """
            (OrderedListAsArray.Cursor) -> None
            Withdraws from the list the object to which this cursor refers.
            """
            if self._offset < 0 \
                    or self._offset >= self._list._count:
                raise IndexError
            if self._list._count == 0:
                raise ContainerEmpty
            i = self._offset
            while i < self._list._count - 1:
                self._list._array[i] = self._list._array[i + 1]
                i += 1
                ++i
            self._list._array[i] = None
            self._list._count -= 1
#}>g

    class Iterator(Iterator):
        """
        Enumerates the items in an ordered list.
        """

        def __init__(self, list):
            """
            (OrderedListAsArray.Iterator, OrderedListAsArray) -> None
            Constructs an iterator for the given ordered list.
            """
            super(OrderedListAsArray.Iterator, self).__init__(list)
            self._position = -1

        def next(self):
            """
            (OrderedListAsArray.Iterator) -> Object
            Returns the next element.
            """
            self._position += 1
            if self._position == self._container.count:
                self._position = -1
                raise StopIteration
            return self._container._array[self._position]

    def __iter__(self):
        """
        (OrderedListAsArray) -> OrderedListAsArray.Iterator
        Returns an iterator for this ordered list.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (OrderedListAsArray, OrderedListAsArray) -> int

        Compares this ordered list with the given ordered list.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "OrderedListAsArray test program."
        print OrderedListAsArray.main.__doc__
        list = OrderedListAsArray(10)
        OrderedList.test(list)
        return 0

if __name__ == "__main__":
    sys.exit(OrderedListAsArray.main(*sys.argv))
