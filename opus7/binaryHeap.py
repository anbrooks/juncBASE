#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: binaryHeap.py,v $
#   $Revision: 1.21 $
#
#   $Id: binaryHeap.py,v 1.21 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BinaryHeap class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.21 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.priorityQueue import PriorityQueue
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.array import Array
from opus7.exception import *

#{
class BinaryHeap(PriorityQueue):
    """
    Binary heap class implemented using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, length = 0):
        """
        (BinaryHeap [,int]) -> None
        Constructs a binary heap of the given length.
        """
        super(BinaryHeap, self).__init__()
        self._array = Array(length, 1) # Base index is 1.

    def purge(self):
        """
        (BinaryHeap) -> None
        Purges this binary heap.
        """
        while self._count > 0:
            self._array[self._count] = None
            self._count -= 1
#}>a

#{
    def enqueue(self, obj):
        """
        (BinaryHeap, Object) -> None
        Enqueues the given object in this binary heap.
        """
        if self._count == len(self._array):
            raise ContainerFull
        self._count += 1
        i = self._count
        while i > 1 and self._array[i/2] > obj:
            self._array[i] = self._array[i / 2]
            i /= 2
        self._array[i] = obj
#}>b

#{
    def getMin(self):
        """
        (BinaryHeap) -> Object
        Returns the object in this binary heap with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._array[1]
#}>c

#{
    def dequeueMin(self):
        """
        (BinaryHeap) -> Object
        Dequeues and returns the object in this binary heap
        with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._array[1]
        last = self._array[self._count]
        self._count -= 1
        i = 1
        while 2 * i < self._count + 1:
            child = 2 * i
            if child + 1 < self._count + 1 \
                and self._array[child + 1] < self._array[child]:
                child += 1
            if last <= self._array[child]:
                break
            self._array[i] = self._array[child]
            i = child
        self._array[i] = last
        return result
#}>d

    def getIsFull(self):
        """
        (BinaryHeap) -> bool
        Returns true if this binary heap is full.
        """
        return self._count == len(self_.array) - 1

    def accept(self, visitor):
        """
        (BinaryHeap, Visitor) -> None
        Makes the given visitor visit all the objects in the binary heap.
        """
        assert isinstance(visitor, Visitor)
        for i in xrange(1, self._count):
            visitor.visit(self._array[i])
            if visitor.isDone:
                return

    class Iterator(Iterator):
        """
        Enumerates the elements of a binary heap.
        """

        def __init__(self, heap):
            """
            (BinaryHeap.Iterator, BinaryHeap) -> None
            Constructs an iterator for the given binary heap.
            """
            super(BinaryHeap.Iterator, self).__init__(heap)
            self._position = 0

        def next(self):
            """
            (BinaryHeap.Iterator) -> Object
            Returns the next object in this binary heap.
            """
            self._position += 1
            if self._position > self.container._count:
                self._position = 0
                raise StopIteration
            return self.container._array[self._position]

    def __iter__(self):
        """
        (BinaryHeap) -> BinaryHeap.Iterator.
        Returns an interator for this binary heap.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (BinaryHeap, BinaryHeap) -> int

        Compares this binary heap with the given binary heap.
        """
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "BinaryHeap test program."
        print BinaryHeap.main.__doc__
        pqueue = BinaryHeap(256)
        PriorityQueue.test(pqueue)
        return 0

if __name__ == "__main__":
    sys.exit(BinaryHeap.main(*sys.argv))
