#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: chainedScatterTable.py,v $
#   $Revision: 1.35 $
#
#   $Id: chainedScatterTable.py,v 1.35 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the ChainedScatterTable class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.35 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object
from opus7.hashTable import HashTable
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class ChainedScatterTable(HashTable):
    """
    Hash table implemented as a chained scatter table using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    NULL = -1

    class Entry(object):
        """
        An entry in a chained scatter table.
        """

        def __init__(self, obj, next):
            """
            (ChainedScatterTable.Entry, Object, int) -> None
            Constructs an entry in a chained scatter table.
            """
            super(ChainedScatterTable.Entry, self).__init__()
            self._obj = obj
            self._next = next
#}>a

#{
    def __init__(self, length):
        """
        (ChainedScatterTable, int) -> None
        Constructs a chained scatter table with the given length.
        """
        super(ChainedScatterTable, self).__init__()
        self._array = Array(length)
        for i in xrange(len(self._array)):
            self._array[i] = self.Entry(None, self.NULL)

    def __len__(self):
        """
        (ChainedScatterTable) -> int
        Returns the length of this chained scatter table.
        """
        return len(self._array)

    def purge(self):
        """
        (ChainedScatterTable) -> None
        Purges this chained scatter table.
        """
        for i in len(self):
            self._array[i] = self.Entry(None, self.NULL)
        self._count = 0
#}>b

    def getIsFull(self):
        """
        (ChainedScatterTable) -> bool
        Returns true if this chained scatter table is full.
        """
        return self._count == len(self)

    def accept(self, visitor):
        """
        (ChainedScatterTable, Visitor) -> None
        Makes the given visitor visit all the objects
        in this chained scatter table.
        """
        assert isinstance(visitor, Visitor)
        for i, entry in enumerate(self._array):
            if entry._obj is not None:
                visitor.visit(entry._obj)
                if visitor.isDone:
                    return

#{
    def insert(self, obj):
        """
        (ChainedScatterTable, Object) -> None
        Inserts the given object into this chained scatter table.
        """
        if self._count == len(self):
            raise ContainerFull
        probe = self.h(obj)
        if self._array[probe]._obj is not None:
            while self._array[probe]._next != self.NULL:
                probe = self._array[probe]._next
            tail = probe
            probe = (probe + 1) % len(self)
            while self._array[probe]._obj is not None:
                probe = (probe + 1) % len(self)
            self._array[tail]._next = probe
        self._array[probe] = self.Entry(obj, self.NULL)
        self._count += 1

    def find(self, obj):
        """
        (ChainedScatterTable, Object) -> Object
        Returns the object in this chained scatter table
       that matches the given object.
       """
        probe = self.h(obj)
        while probe != self.NULL:
            if obj == self._array[probe]._obj:
                return self._array[probe]._obj
            probe = self._array[probe]._next
        return None
#}>c

#{
    def withdraw(self, obj):
        """
        (ChainedScatterTable, Object) -> None
        Withdraws the given object from this chained scatter table.
        """
        if self._count == 0:
            raise ContainerEmpty
        i = self.h(obj)
        while i != self.NULL and self._array[i]._obj is not obj:
            i = self._array[i]._next
        if i == self.NULL:
            raise KeyError
        while True:
            j = self._array[i]._next
            while j != self.NULL:
                h = self.h(self._array[j]._obj)
                contained = False
                k = self._array[i]._next
                while k != self._array[j]._next \
                        and not contained:
                    if k == h:
                        contained = True
                    k = self._array[k]._next
                if not contained:
                    break
                j = self._array[j]._next
            if j == self.NULL:
                break
            self._array[i]._obj = self._array[j]._obj
            i = j
        self._array[i] = self.Entry(None, self.NULL)
        j = (i + len(self) - 1) % len(self)
        while j != i:
            if self._array[j]._next == i:
                self._array[j]._next = self.NULL
                break
            j = (j + len(self) - 1) % len(self)
        self._count -= 1
#}>d

    def __contains__(self, obj):
        """
        (ChainedScatterTable, Object) -> bool
        Returns true if the given object is in this chained scatter table.
        """
        probe = self.h(obj)
        while probe != self.NULL:
            if array[probe]._obj is obj:
                return True
            probe = array[probe]._next
        return False

    class Iterator(Iterator):
        """
        Enumerates the elements of a chained scatter table.
        """

        def __init__(self, table):
            """
            (ChainedScatterTable.Iterator, ChainedScatterTable) -> None
            Constructs an interator for the given chained scatter table.
            """
            super(ChainedScatterTable.Iterator, self).__init__(table)
            self._position = -1

        def next(self):
            """
            (ChainedScatterTable.Iterator) -> Object
            Returns the next object in the chained scatter table.
            """
            self._position += 1
            while self._position < len(self._container._array):
                if self._container._array[self._position]._obj is not None:
                    break
                self._position += 1
            if self._position == len(self._container._array):
                self._position = -1
                raise StopIteration
            return self._container._array[self._position]._obj


    def __iter__(self):
        """
        (ChainedScatterTable) -> ChainedScatterTable.Iterator
        Returns an iterator for this chained scatter table.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (ChainedScatterTable, ChainedScatterTable) -> int

        Compares this scatter table with the given scatter table.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "ChainedScatterTable test program."
        print ChainedScatterTable.main.__doc__
        hashTable = ChainedScatterTable(57)
        HashTable.test(hashTable)
        return 0

if __name__ == "__main__":
    sys.exit(ChainedScatterTable.main(*sys.argv))
