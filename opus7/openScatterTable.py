#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: openScatterTable.py,v $
#   $Revision: 1.29 $
#
#   $Id: openScatterTable.py,v 1.29 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the OpenScatterTable class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.29 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object
from opus7.hashTable import HashTable
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class OpenScatterTable(HashTable):
    """
    Hash table implemented as an open scatter table using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    EMPTY = 0
    OCCUPIED = 1
    DELETED = 2

    class Entry(object):
        """
        An entry in an open scatter table.
        """

        def __init__(self, state, obj):
            """
            (OpenScatterTable.Entry, int, Object) -> None
            Constructs an entry in an open scatter table.
            """
            super(OpenScatterTable.Entry, self).__init__()
            self._state = state
            self._obj = obj
#}>a

#{
    def __init__(self, length):
        """
        (OpenScatterTable, int) -> None
        Constructs an open scatter table of the given length.
        """
        super(OpenScatterTable, self).__init__()
        self._array = Array(length)
        for i in xrange(len(self._array)):
            self._array[i] = self.Entry(self.EMPTY, None)

    def __len__(self):
        """
        (OpenScatterTable) -> None
        Returns the length of this open scatter table.
        """
        return len(self._array)

    def purge(self):
        """
        (OpenScatterTable) -> None
        Purges this open scatter table.
        """
        for i in xrange(len(self._array)):
            self._array[i] = self.Entry(self.EMPTY, None)
        self._count = 0
#}>b

    def getIsFull(self):
        """
        (OpenScatterTable) -> bool
        Returns true if this open scatter table is full.
        """
        return self._count == len(self)

    def accept(self, visitor):
        """
        (OpenScatterTable, Visitor) -> None
        Makes the given visitor visit all the objects
        in this open scatter table.
        """
        assert isinstance(visitor, Visitor)
        for i, entry in enumerate(self._array):
            if entry._state == self.OCCUPIED:
                visitor.visit(entry._obj)
                if visitor.isDone:
                    return

#{
    def c(self, i):
        """
        (OpenScatterTable, int) -> int
        Linear probing function.
        """
        return i

    def findUnoccupied(self, obj):
        """
        (OpenScatterTable, Object) -> int
        Returns the index of an unoccupied entry in this open scatter table.
        """
        hash = self.h(obj)
        for i in xrange(self._count + 1):
            probe = (hash + self.c(i)) % len(self)
            if self._array[probe]._state != self.OCCUPIED:
                return probe
        raise ContainerFull

    def insert(self, obj):
        """
        (OpenScatterTable, Object) -> None
        Inserts the given object into this open scatter table.
        """
        if self._count == len(self):
            raise ContainerFull
        offset = self.findUnoccupied(obj)
        self._array[offset] = self.Entry(self.OCCUPIED, obj)
        self._count += 1
#}>c

#{
    def findMatch(self, obj):
        """
        (OpenScatterTable, Object) -> int
        Finds the index of an object in this open scatter table
       that matches the given object.
       """
        hash = self.h(obj)
        for i in xrange(len(self._array)):
            probe = (hash + self.c(i)) % len(self)
            if self._array[probe]._state == self.EMPTY:
                break
            if self._array[probe]._state == self.OCCUPIED and \
                    self._array[probe]._obj == obj:
                return probe
        return -1

    def find(self, obj):
        """
        (OpenScatterTable, Object) -> Object
        Returns the object in this open scatter table
       that matches the given object.
       """
        offset = self.findMatch(obj)
        if offset >= 0:
            return self._array[offset]._obj
        else:
            return None
#}>d

    def findInstance(self, obj):
        """
        (OpenScatterTable, Object) -> int
        Finds the index of the given object in this open scatter table.
        """
        hash = self.h(obj)
        for i in xrange(len(self._array)):
            probe = (hash + self.c(i)) % len(self)
            if self._array[probe]._state == self.EMPTY:
                break
            if self._array[probe]._state == self.OCCUPIED and \
                self._array[probe]._obj is obj:
                return probe
        return -1

#{
    def withdraw(self, obj):
        """
        (OpenScatterTable, Object) -> None
        Withdraws the given object from this open scatter table.
        """
        if self._count == 0:
            raise ContainerEmpty
        offset = self.findInstance(obj)
        if offset < 0:
            raise KeyError
        self._array[offset] = self.Entry(self.DELETED, None)
        self._count -= 1
#}>e

    def __contains__(self, obj):
        """
        (OpenScatterTable, Object) -> bool
        Returns true if the given object is in this open scatter table.
        """
        return self.findInstance(obj) >= 0

    class Iterator(Iterator):
        """
        Enumerates the elements of an open scatter table.
        """

        def __init__(self, table):
            """
            (OpenScatterTable.Iterator, OpenScatterTable) -> None
            Constructs an interator for the given open scatter table.
            """
            super(OpenScatterTable.Iterator, self).__init__(table)
            self._position = -1

        def next(self):
            """
            (OpenScatterTable.Iterator) -> Object
            Returns the next object in the open scatter table.
            """
            self._position += 1
            while self._position < len(self._container):
                if self._container._array[self._position]._state == \
                        self._container.OCCUPIED:
                    break
                self._position += 1
            if self._position == len(self._container):
                self._position = -1
                raise StopIteration
            return self._container._array[self._position]._obj
        
    def __iter__(self):
        """
        (OpenScatterTable) -> OpenScatterTable.Iterator
        Returns an interator for this open scatter table.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (OpenScatterTable, OpenScatterTable) -> int

        Compares this scatter table with the given scatter table.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "OpenScatterTable test program."
        print OpenScatterTable.main.__doc__
        hashTable = OpenScatterTable(57)
        HashTable.test(hashTable)
        return 0

if __name__ == "__main__":
    sys.exit(OpenScatterTable.main(*sys.argv))
