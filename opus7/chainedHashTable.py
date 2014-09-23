#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: chainedHashTable.py,v $
#   $Revision: 1.29 $
#
#   $Id: chainedHashTable.py,v 1.29 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the ChainedHashTable class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.29 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.hashTable import HashTable
from opus7.array import Array
from opus7.linkedList import LinkedList
from opus7.iterator import Iterator
from opus7.visitor import Visitor

#{
class ChainedHashTable(HashTable):
    """
    Hash table implemented using an array of linked lists.
    """

#}@head

#{
    # ...
#}@tail

#{
    def __init__(self, length):
        """
        (ChainedHashTable, int) -> None
        Constructs a chained hash table with the given length.
        """
        super(ChainedHashTable, self).__init__()
        self._array = Array(length)
        for i in xrange(len(self._array)):
            self._array[i] = LinkedList()

    def __len__(self):
        """
        (ChainedHashTable) -> int
        Returns the length of this chained hash table.
        """
        return len(self._array)

    def purge(self):
        """
        (ChainedHashTable) -> None
        Purges this chained hash table.
        """
        for i in xrange(len(self._array)):
            self._array[i].purge()
        self._count = 0
#}>a

#{
    def insert(self, obj):
        """
        (ChainedHashTable, Object) -> None
        Inserts the given object into this chained hash table.
        """
        self._array[self.h(obj)].append(obj)
        self._count += 1

    def withdraw(self, obj):
        """
        (ChainedHashTable, Object) -> None
        Withdraws the given object from this chained hash table.
        """
        self._array[self.h(obj)].extract(obj)
        self._count -= 1
#}>b

    def __contains__(self, obj):
        """
        (ChainedHashTable, Object) -> bool
        Returns true if the given object is in this chained hash table.
        """
        ptr = self._array[self.h(obj)].head
        while ptr is not None:
            if ptr.datum is obj:
                return True
            ptr = ptr.next
        return False

#{
    def find(self, obj):
        """
        (ChainedHashTable, Object) -> Object
        Returns the obj in this hash table
       that matches the given object.
       """
        ptr = self._array[self.h(obj)].head
        while ptr is not None:
            datum = ptr.datum
            if obj == datum:
                return datum
            ptr = ptr.next
        return None
#}>c

    def accept(self, visitor):
        """
        (ChainedHashTable, Visitor) -> None
        Makes the given visitor visit all the objects in this hash table.
        """
        assert isinstance(visitor, Visitor)
        for i, entry in enumerate(self._array):
            ptr = entry.head
            while ptr is not None:
                visitor.visit(ptr.datum)
                if visitor.isDone:
                    return
                ptr = ptr.next

    class Iterator(Iterator):
        """
        Enumerates the._elements of a chained hash table.
        """

        def __init__(self, hashTable):
            """
            (ChainedHashTable.Iterator, ChainedHashTable) -> None
            Constructs an iterator for the given hash table.
            """
            super(ChainedHashTable.Iterator, self).__init__(hashTable)
            self._element = None
            self._position = -1

        def next(self):
            """
            (ChainedHashTable.Iterator) -> Object
            Returns the next object in the hash table.
            """
            if self._element is not None:
                self._element = self._element.next
            if self._element is None:
                self._position += 1
                while self._position < len(self._container):
                    self._element = self._container._array[self._position].head
                    if self._element is not None:
                        break
                    self._position += 1
                if self._position == len(self._container):
                    self._position = -1
            if self._element is None:
                raise StopIteration
            return self._element.datum

    def __iter__(self):
        """
        (ChainedHashTable) -> ChainedHashTable.Iterator
        Returns an iterator for this hash table.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (ChainedHashTable, ChainedHashTable) -> int

        Compares this hash table with the given hash table.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "ChainedHashTable test program."
        print ChainedHashTable.main.__doc__
        hashTable = ChainedHashTable(57)
        HashTable.test(hashTable)
        return 0

if __name__ == "__main__":
    sys.exit(ChainedHashTable.main(*sys.argv))
