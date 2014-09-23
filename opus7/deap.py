#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: deap.py,v $
#   $Revision: 1.22 $
#
#   $Id: deap.py,v 1.22 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Deap class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.22 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.doubleEndedPriorityQueue import DoubleEndedPriorityQueue
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.array import Array
from opus7.exception import *

class Deap(DoubleEndedPriorityQueue):
    """
    Double-ended heap.
    """

    def __init__(self, length = 0):
        """
        (Deap, int) -> None
        Constructs a deap with the given length.
        """
        super(Deap, self).__init__()
        self._array = Array(length + 1, 1)

    def purge(self):
        """
        (Deap) -> None
        Purges this deap.
        """
        self._array = Array(len(self_.array), 1)
        self._count = 0

    @staticmethod
    def log2(i):
        """
        (int) -> int
        Returns ceil(log_2(i)).
        """
        result = 0
        while (1 << result) <= i:
            result += 1
        return result - 1

    @staticmethod
    def mask(i):
        """
        (int) -> int
        Returns 2^(ceil(log2(i)-1).
        """
        return 1 << (Deap.log2(i) - 1)

    def dual(self, i):
        """
        (Deap, i) -> int
        Returns the position of the dual of the given position.
        """
        m = Deap.mask(i)
        result = i ^ m
        if (result & m) != 0:
            if result >= self._count + 2:
                result /= 2
        else:
            if 2 * result < self._count + 2:
                result *= 2
                if result + 1 < self._count + 2 \
                        and self._array[result + 1] > self._array[result]:
                    result += 1
        return result

    def getMin(self):
        """
        (Deap) -> Object
        Returns the object in this deap with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._array[2]

    def getMax(self):
        """
        (Deap) -> Object
        Returns the object in this deap with the largest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        if self._count == 1:
            return self._array[2]
        else:
            return self._array[3]

    def insertMin(self, pos, obj):
        """
        (Deap, int, Object) -> None
        Inserts the given object into the min heap of this deap
        starting from the given position.
        """
        i = pos
        while i > 2 and self._array[i / 2] > obj:
            self._array[i] = self._array[i / 2]
            i /= 2
        self._array[i] = obj

    def insertMax(self, pos, obj):
        """
        (Deap, int, Object) -> None
        Inserts the given object into the max heap of this deap
        starting from the given position.
        """
        i = pos
        while i > 3 and self._array[i / 2] < obj:
            self._array[i] = self._array[i / 2]
            i /= 2
        self._array[i] = obj

    def enqueue(self, obj):
        """
        (Deap, Object) -> None
        Enqueues the given object in this deap.
        """
        if self._count == len(self._array) - 1:
            raise ContainerFull
        self._count += 1
        if self._count == 1:
            self._array[2] = obj
        else:
            i = self._count + 1
            j = self.dual(i)
            if (i & Deap.mask(i)) != 0:
                if obj >= self._array[j]:
                    self.insertMax(i, obj)
                else:
                    self._array[i] = self._array[j]
                    self.insertMin(j, obj)
            else:
                if obj < self._array[j]:
                    self.insertMin(i, obj)
                else:
                    self._array[i] = self._array[j]
                    self.insertMax(j, obj)

    def dequeueMin(self):
        """
        (Deap) -> Object
        Dequeues and returns the object in this deap
        with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._array[2]
        last = self._array[self._count + 1]
        self._count -= 1
        if self._count <= 1:
            self._array[2] = last
        else:
            i = 2
            while 2 * i < self._count + 2:
                child = 2 * i
                if child + 1 < self._count + 2 \
                        and self._array[child + 1] < self._array[child]:
                    child += 1
                self._array[i] = self._array[child]
                i = child
            j = self.dual(i)
            if last <= self._array[j]:
                self.insertMin(i, last)
            else:
                self._array[i] = self._array[j]
                self.insertMax(j, last)
        return result

    def dequeueMax(self):
        """
        (Deap) -> Object
        Dequeues and returns the object in this deap
        with the largest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        if self._count == 1:
            self._count -= 1
            return self._array[2]
        elif self._count == 2:
            self._count -= 1
            return self._array[3]
        else:
            result = self._array[3]
            last = self._array[self._count + 1]
            self._count -= 1
            i = 3
            while 2 * i < self._count + 2:
                child = 2 * i
                if child + 1 < self._count + 2 \
                        and self._array[child + 1] > self._array[child]:
                    child += 1
                self._array[i] = self._array[child]
                i = child
            j = self.dual(i)
            if last >= self._array[j]:
                self.insertMax(i, last)
            else:
                self._array[i] = self._array[j]
                self.insertMin(j, last)
            return result

    def getIsFull(self):
        """
        (Deap) -> bool
        Returns true if this deap is full
        """
        return self._count == len(self._array) - 2

    def accept(self, visitor):
        """
        (Deap, Visitor) -> None
        Makes the given visitor visit all the objects in this deap.
        """
        assert isinstance(visitor, Visitor)
        for i in xrange(2, self._count + 2):
            visitor.visit(self._array[i])
            if visitor.isDone:
                return

    class Iterator(Iterator):
        """
        Enumerates the objects in a deap.
        """

        def __init__(self, deap):
            """
            (Deap.Iterator, Deap) -> None
            Constructs an iterator for the given deap.
            """
            super(Deap.Iterator, self).__init__(deap)
            self._position = 1

        def next(self):
            """
            (Deap.Iterator) -> Object
            Returns the next object in the deap.
            """
            self._position += 1
            if self._position == self._container._count - 2:
                self._position = 1
                raise StopIteration
            return self._container._array[self._position]


    def __iter__(self):
        """
        (Deap) -> Deap.Iterator
        Returns an iterator for this deap.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (Deap, Deap) -> int

        Compares this deap with the given deap.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "Deap test program."
        print Deap.main.__doc__
        depqueue = Deap(256)
        DoubleEndedPriorityQueue.test(depqueue)
        return 0

if __name__ == "__main__":
    sys.exit(Deap.main(*sys.argv))
