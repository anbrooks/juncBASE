#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: setAsBitVector.py,v $
#   $Revision: 1.30 $
#
#   $Id: setAsBitVector.py,v 1.30 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SetAsBitVector class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.30 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import warnings
from opus7.set import Set
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor

# Suppress warnings about left shifts.
warnings.filterwarnings("ignore", "x<<y", FutureWarning)

#{
class SetAsBitVector(Set):
    """
    Set implemented using a bit vector.
    """

#}@head

#{

    # ...
#}@tail

#{
    BITS = 32 # Number of bits in an integer.

    def __init__(self, n):
        """
        (SetAsBitVector, int) -> None
        Constructs a set with the given universe size.
        """
        super(SetAsBitVector, self).__init__(n)
        self._vector = Array((n + self.BITS - 1) / self.BITS)
        for i in xrange(len(self._vector)):
            self._vector[i] = 0
#}>a

#{
    def insert(self, item):
        """
        (SetAsBitVector, int) -> None
        Inserts the given element into this set.
        """
        self._vector[item / self.BITS] |= 1 << item % self.BITS

    def withdraw(self, item):
        """
        (SetAsBitVector, int) -> None
        Withdraws the given element from this set.
        """
        self._vector[item / self.BITS] &= \
            ~(1 << item % self.BITS)

    def __contains__(self, item):
        """
        (SetAsBitVector, int) -> bool
        Returns true if the given element is in this set.
        """
        return (self._vector[item / self.BITS] \
            & (1 << item % self.BITS)) != 0
#}>b

    def getIsEmpty(self):
        """
        (SetAsBitVector) -> bool
        Returns true if this set is empty.
        """
        for item in xrange(self._universeSize):
            if item in self:
                return False
        return True

    def getIsFull(self):
        """
        (SetAsBitVector) -> bool
        Returns true if this set is full.
        """
        for item in xrange(self._universeSize):
            if item not in self:
                return False
        return True

    def purge(self):
        """
        (SetAsBitVector) -> None
        Purges this set.
        """
        for i in xrange(len(self._vector)):
            self._vector[i] = 0

    def accept(self, visitor):
        """
        (SetAsBitVector, Visitor) -> None
        Makes the given visitor visit all the elements in this set.
        """
        assert isinstance(visitor, Visitor)
        for item in xrange (self._universeSize):
            if item in self:
                visitor.visit(item)

    def getCount(self):
        """
        (SetAsBitVector) -> int
        Returns the number of elements in this set.
        """
        result = 0
        for item in xrange(self._universeSize):
            if item in self:
                result += 1
        return result

#{
    def __or__(self, set):
        """
        (SetAsBitVector, SetAsBitVector) -> SetAsBitVector
        Returns the union of this set an the given set.
        """
	assert isinstance(set, SetAsBitVector)
	assert self._universeSize == set._universeSize
        result = SetAsBitVector(self._universeSize)
        for i in xrange(len(self._vector)):
            result._vector[i] = self._vector[i] | set._vector[i]
        return result

    def __and__(self, set):
        """
        (SetAsBitVector, SetAsBitVector) -> SetAsBitVector
        Returns the intersection of this set an the given set.
        """
	assert isinstance(set, SetAsBitVector)
	assert self._universeSize == set._universeSize
        result = SetAsBitVector(self._universeSize)
        for i in xrange(len(self._vector)):
            result._vector[i] = self._vector[i] & set._vector[i]
        return result

    def __sub__(self, set):
        """
        (SetAsBitVector, SetAsBitVector) -> SetAsBitVector
        Returns the difference of this set an the given set.
        """
	assert isinstance(set, SetAsBitVector)
	assert self._universeSize == set._universeSize
        result = SetAsBitVector(self._universeSize)
        for i in xrange(len(self._vector)):
            result._vector[i] = self._vector[i] & ~set._vector[i]
        return result
#}>c

    def __le__(self, set):
        """
        (SetAsBitVector, SetAsBitVector) -> bool
        Returns true if this set is a subset of the given set.
        """
	assert isinstance(set, SetAsBitVector)
	assert self._universeSize == set._universeSize
        for i in xrange(len(self._vector)):
            if (self._vector[i] & ~set._vector[i]) != 0:
                return False
        return True

    def __eq__(self, set):
        """
        (SetAsBitVector, SetAsBitVector) -> bool
        Returns true if this set equals the given set.
        """
	assert isinstance(set, SetAsBitVector)
	assert self._universeSize == set._universeSize
        for i in xrange(len(self._vector)):
            if self._vector[i] != set._vector[i]:
                return False
        return True

    class Iterator(Iterator):
        """
        Enumerates the elements of a SetAsBitVector.
        """

        def __init__(self, set):
            """
            (SetAsBitVector.Iterator, SetAsBitVector) -> None
            Constructs an interator for the given set.
            """
            super(SetAsBitVector.Iterator, self).__init__(set)
            self._item = -1
        
        def next(self):
            """
            (SetAsBitVector.Iterator) -> Object
            Returns the next element in the set.
            """
            self._item += 1
            while self._item < self._container._universeSize:
                if self._item in self._container:
                    break
                self._item += 1
            if self._item == self._container._universeSize:
                self._item = -1
                raise StopIteration
            return self._item

    def __iter__(self):
        """
        (SetAsBitVector) -> SetAsBitVector.Iterator
        Returns an interator for this set.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (SetAsBitVector, SetAsBitVector) -> int

        Compares this set with the given set.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "SetAsBitVector test program."
        print SetAsBitVector.main.__doc__
        Set.test(SetAsBitVector (57), SetAsBitVector(57), SetAsBitVector(57))
        return 0

if __name__ == "__main__":
    sys.exit(SetAsBitVector.main(*sys.argv))
