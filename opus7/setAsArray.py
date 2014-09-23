#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: setAsArray.py,v $
#   $Revision: 1.28 $
#
#   $Id: setAsArray.py,v 1.28 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SetAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.28 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.set import Set
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor

#{
class SetAsArray(Set):
    """
    Set implemented using an array of boolean values.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, n):
        """
        (SetAsArray, int) -> None
        Constructs a set with the given universe size.
        """
        super(SetAsArray, self).__init__(n)
        self._array = Array(self._universeSize)
        for item in xrange(0, self._universeSize):
            self._array[item] = False
#}>a

#{
    def insert(self, item):
        """
        (SetAsArray, int) -> None
        Inserts the given element into this set.
        """
        self._array[item] = True

    def __contains__(self, item):
        """
        (SetAsArray, int) -> bool
        Returns true if the given element is in this set.
        """
        return self._array[item]

    def withdraw(self, item):
        """
        (SetAsArray, int) -> None
        Withdraws the given element from this set.
        """
        self._array[item] = False
#}>b

    def getIsEmpty(self):
        """
        (SetAsArray) -> bool
        Returns true if this set is empty.
        """
        for item in xrange(0, self._universeSize):
            if item in self:
                return False
        return True

    def getIsFull(self):
        """
        (SetAsArray) -> bool
        Returns true if this set is full.
        """
        for item in xrange(0, self._universeSize):
            if item not in self:
                return False
        return True

    def purge(self):
        """
        (SetAsArray) -> None
        Purges this set.
        """
        for item in xrange(0, self._universeSize):
            self._array[item] = False

    def accept(self, visitor):
        """
        (SetAsArray, Visitor) -> None
        Makes the given visitor visit all the elements of this set.
        """
        assert isinstance(visitor, Visitor)
        for item in xrange(0, self._universeSize):
            if item in self:
                visitor.visit(item)

    def getCount(self):
        """
        (SetAsArray) -> None
        Returns the number of elements in this set.
        """
        result = 0
        for item in xrange(0, self._universeSize):
            if item in self:
                result += 1
        return result

#{
    def __or__(self, set):
        """
        (SetAsArray, SetAsArray) -> None
        Returns the union of this set and the given set.
        """
	assert isinstance(set, SetAsArray)
	assert self._universeSize == set.universeSize
        result = SetAsArray(self._universeSize)
        for i in xrange(0, self._universeSize):
            result._array[i] = self._array[i] or set._array[i]
        return result

    def __and__(self, set):
        """
        (SetAsArray, SetAsArray) -> None
        Returns the intersection of this set and the given set.
        """
	assert isinstance(set, SetAsArray)
	assert self._universeSize == set.universeSize
        result = SetAsArray(self._universeSize)
        for i in xrange(0, self._universeSize):
            result._array[i] = self._array[i] and set._array[i]
        return result

    def __sub__(self, set):
        """
        (SetAsArray, SetAsArray) -> None
        Returns the difference of this set and the given set.
        """
	assert isinstance(set, SetAsArray)
	assert self._universeSize == set.universeSize
        result = SetAsArray(self._universeSize)
        for i in xrange(0, self._universeSize):
            result._array[i] = self._array[i] \
                and not set._array[i]
        return result
#}>c

#{
    def __eq__(self, set):
        """
        (SetAsArray, SetAsArray) -> bool
        Returns true if this set is equal to the given set.
        """
	assert isinstance(set, SetAsArray)
	assert self._universeSize == set.universeSize
        for i in xrange(0, self._universeSize):
            if self._array[i] != set._array[i]:
                return False
        return True

    def __le__(self, set):
        """
        (SetAsArray, SetAsArray) -> bool
        Returns true if this set is a proper subset to the given set.
        """
	assert isinstance(set, SetAsArray)
	assert self._universeSize == set.universeSize
        for i in xrange(0, self._universeSize):
            if self_array[i] and not set._array[i]:
                return False
        return True
#}>d

    class Iterator(Iterator):
        """
        Enumerates the elements of a SetAsArray.
        """

        def __init__(self, set):
            """
            (SetAsArray.Iterator, SetAsArray) -> None
            Constructs an enumerator for the given set.
            """
            super(SetAsArray.Iterator, self).__init__(set)
            self._item = -1
        
        def next(self):
            """
            (SetAsArray.Iterator) -> int
            Returns the next element in the set.
            """
            self._item += 1
            while self._item < self._container.universeSize:
                if self._item in self._container:
                    break
                self._item += 1
            if self._item == self._container.universeSize:
                self._item = -1
                raise StopIteration
            return self._item
        
    def __iter__(self):
        """
        (SetAsArray) -> SetAsArray.Iterator
        Returns an interator for this set.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (SetAsArray, SetAsArray) -> int

        Compares this set with the given set.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "SetAsArray test program."
        print SetAsArray.main.__doc__
        Set.test(SetAsArray(32), SetAsArray(32), SetAsArray(32))
        return 0

if __name__ == "__main__":
    sys.exit(SetAsArray.main(*sys.argv))
