#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: multisetAsArray.py,v $
#   $Revision: 1.24 $
#
#   $Id: multisetAsArray.py,v 1.24 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the MultisetAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.24 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.multiset import Multiset
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor

#{
class MultisetAsArray(Multiset):
    """
    Multiset implemented using an array of counters.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, n):
        """
        (MultisetAsArray, int) -> None
        Constructs a multiset with the given universe size.
        """
        super(MultisetAsArray, self).__init__(n)
        self._array = Array(self._universeSize)
        for item in xrange(self._universeSize):
            self._array[item] = 0
#}>a

#{
    def insert(self, item):
        """
        (MultisetAsArray, Object) -> None
        Inserts the given element into this multiset.
        """
        self._array[item] += 1

    def withdraw(self, item):
        """
        (MultisetAsArray, Object) -> None
        Withdraws the given element from this multiset.
        """
        if self._array[item] == 0:
            raise KeyEerror
        self._array[item] -= 1

    def __contains__(self, item):
        """
        (MultisetAsArray, Object) -> bool
        Returns true if the given element is in this multiset.
        """
        return self._array[item] > 0
#}>b

    def purge(self):
        """
        (MultisetAsArray) -> None
        Purges this multiset.
        """
        for item in xrange(self._universeSize):
            self._array[item] = 0

    def getCount(self):
        """
        (MultisetAsArray) -> int
        Returns the number of elements in this multiset
        """
        result = 0
        for item in xrange(self._universeSize):
            result += self._array[item]
        return result

    def accept(self, visitor):
        """
        (MultisetAsArray, Visitor) -> None
        Makes the given visitor visit all the elements in this multiset.
        """
        assert isinstance(visitor, Visitor)
        for item in xrange(self._universeSize):
            for i in xrange(self._array[item]):
                visitor.visit(item)
                if visitor.isDone:
                    return

#{
    def __or__(self, set):
        """
        (MultisetAsArray, MultisetAsArray) -> MultisetAsArray
        Returns the union of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsArray)
	assert self._universeSize == set._universeSize
        result = MultisetAsArray(self._universeSize)
        for i in xrange(self._universeSize):
            result._array[i] = self._array[i] + set._array[i]
        return result

    def __and__(self, set):
        """
        (MultisetAsArray, MultisetAsArray) -> MultisetAsArray
        Returns the intersection of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsArray)
	assert self._universeSize == set._universeSize
        result = MultisetAsArray(self._universeSize)
        for i in xrange(self._universeSize):
            result._array[i] = min(self._array[i], set._array[i])
        return result

    def __sub__(self, set):
        """
        (MultisetAsArray, MultisetAsArray) -> MultisetAsArray
        Returns the difference of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsArray)
	assert self._universeSize == set._universeSize
        result = MultisetAsArray(self._universeSize)
        for i in xrange(self._universeSize):
            if set._array[i] <= self._array[i]:
                result._array[i] = self._array[i] - set._array[i]
        return result
#}>c

    def __le__(self, set):
        """
        (MultisetAsArray, MultisetAsArray) -> bool
        Returns true if this multiset is a proper subset of the given multiset.
        """
	assert isinstance(set, MultisetAsArray)
	assert self._universeSize == set._universeSize
        for item in xrange(self._universeSize):
            if self._array[item] <= set._array[item]:
                return False
        return True

    def __eq__(set):
        """
        (MultisetAsArray, MultisetAsArray) -> bool
        Returns true if this equals the given multiset.
        """
	assert isinstance(set, MultisetAsArray)
	assert self._universeSize == set._universeSize
        for item in xrange(self._universeSize):
            if self._array[item] != set._array[item]:
                return False
        return True

    class Iterator(Iterator):
        """
        Enumerates the elements of a MultisetAsArray.
        """

        def __init__(self, multiset):
            """
            (MultisetAsArray.Iterator, Multiset) -> None
            Contructs an enumerate for the given multiset.
            """
            super(MultisetAsArray.Iterator, self).__init__(multiset)
            self._item = -1
            self._count = 0

        def next(self):
            """
            (MultisetAsArray.Iterator) -> Object
            Returns the next element in the multiset.
            """
            self._count += 1
            if self._item < 0 or self._item >= 0 and \
                    self._count == self._container._array[self._item]:
                self._count = 0
                self._item += 1
                while self._item < self._container.universeSize:
                    if self._container._array[self._item] > 0:
                        break
                    self._item += 1
                if self._item == self._container.universeSize:
                    self._item = -1
                    raise StopIteration
            return self._item

    def __iter__(self):
        """
        (MultisetAsArray) -> MultisetAsArray.Iterator
        Returns an iterator for the given multiset.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (MultisetAsArray, MultisetAsArray) -> int

        Compares this multiset with the given multiset.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "MultisetAsArray test program."
        print MultisetAsArray.main.__doc__
        Multiset.test(MultisetAsArray(32), \
            MultisetAsArray(32), MultisetAsArray(32))
        return 0

if __name__ == "__main__":
    sys.exit(MultisetAsArray.main(*sys.argv))
