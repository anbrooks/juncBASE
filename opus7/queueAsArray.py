#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: queueAsArray.py,v $
#   $Revision: 1.29 $
#
#   $Id: queueAsArray.py,v 1.29 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the QueueAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.29 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.queue import Queue
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class QueueAsArray(Queue):
    """
    Queue implemented using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size = 0):
        """
        (QueueAsArray [, int]) -> None
        Constructs a queue of the given size.
        """
        super(QueueAsArray, self).__init__()
        self._array = Array(size)
        self._head = 0
        self._tail = size - 1

    def purge(self):
        """
        (QueueAsArray) -> None
        Purges this queue.
        """
        while self._count > 0:
            self._array[self._head] = None
            self._head = self._head + 1
            if self._head == len(self._array):
                self._head = 0
            self._count -= 1
#}>a

#{
    def getHead(self):
        """
        (QueueAsArray) -> Object
        Returns the object at the head of this queue.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._array[self._head]

    def enqueue(self, obj):
        """
        (QueueAsArray, Object) -> None
        Enqueues the given object to the tail of this queue.
        """
        if self._count == len(self._array):
            raise ContainerFull
        self._tail = self._tail + 1
        if self._tail == len(self._array):
            self._tail = 0
        self._array[self._tail] = obj
        self._count += 1

    def dequeue(self):
        """
        (QueueAsArray) -> Object
        Dequeues the object at the head of this queue.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._array[self._head]
        self._array[self._head] = None
        self._head = self._head + 1
        if self._head == len(self._array):
            self._head = 0
        self._count -= 1
        return result
#}>b

    def getIsFull(self):
        """
        (QueueAsArray) -> bool
        Returns true of this queue is full.
        """
        return self._count == len(self._array)

    def accept(self, visitor):
        """
        (QueueAsArray, Visitor) -> None
        Makes the given visitor visit all the objects in this queue.
        """
        assert isinstance(visitor, Visitor)
        pos = self._head
        for i in xrange(self._count):
            visitor.visit(self._array[pos])
            if visitor.isDone:
                return
            pos = pos + 1
            if ++pos == len(self._array):
                pos = 0

    class Iterator(Iterator):
        "Enumerates the elements of a QueueAsArray."

        def __init__(self, queue):
            """
            (QueueAsArray.Iterator, QueueAsArray) -> None
            Constructs an iterator for the given queue.
            """
            super(QueueAsArray.Iterator).__init__(queue)
            self._position = -1

        def next(self):
            """
            (QueueAsArray.Iterator) -> Object
            Returns the next element.
            """
            self._position = self._position + 1
            if self._position == self._container.getCount():
                self._position = -1
                raise StopIteration
            return queue._array[
                (self._container._head + self._position)
                % len(self._container._array)]

    def __iter__(self):
        """
        (QueueAsArray) -> QueueAsArray.Iterator
        Returns an iterator for this queue.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (QueueAsArray, QueueAsArray) -> int

        Comparse this queue with the given queue.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "QueueAsArray test program."
        print QueueAsArray.main.__doc__
        queue1 = QueueAsArray(5)
        Queue.test(queue1)
        return 0

if __name__ == "__main__":
    sys.exit(QueueAsArray.main(*sys.argv))
