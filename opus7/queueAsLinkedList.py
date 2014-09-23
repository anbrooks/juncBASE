#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: queueAsLinkedList.py,v $
#   $Revision: 1.28 $
#
#   $Id: queueAsLinkedList.py,v 1.28 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the QueueAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.28 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.queue import Queue
from opus7.linkedList import LinkedList
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class QueueAsLinkedList(Queue):
    """
    Queue implemented using a linked list.
    """

#}@head

#{

    # ...
#}@tail


#{
    def __init__(self):
        """
        (QueueAsLinkedList) -> None
        Constructs a queue.
        """
        super(QueueAsLinkedList, self).__init__()
        self._list = LinkedList()

    def purge(self):
        """
        (QueueAsLinkedList) -> None
        Purges this queue.
        """
        self._list.purge()
        self._count = 0
#}>a

#{
    def getHead(self):
        """
        (QueueAsLinkedList) -> Object
        Returns the object at the head of this queue.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._list.first

    def enqueue(self, obj):
        """
        (QueueAsLinkedList, Object) -> None
        Enqueues the given object to the tail of this queue.
        """
        self._list.append(obj)
        self._count += 1

    def dequeue(self):
        """
        (QueueAsLinkedList) -> Object
        Dequeues the object at the head of this queue.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._list.first
        self._list.extract(result)
        self._count -= 1
        return result
#}>b

    def accept(self, visitor):
        """
        (QueueAsLinkedList, Visitor) -> None
        Makes the given visitor visit all the objects in this queue.
        """
        assert isinstance(visitor, Visitor)
        ptr = self._list.head
        while ptr is not None:
            visitor.visit(ptr.datum)
            if visitor.isDone:
                return
            ptr = ptr.next

    class Iterator(Iterator):
        """
        Enumerates the elements of a QueueAsLinkedList.
        """

        def __init__(self, queue):
            """
            (QueueAsLinkedList.Iterator, QueueAsLinkedList) -> None
            Constructs an iterator for the given queue.
            """
            super(QueueAsLinkedList.Iterator).__init__(queue)
            self._position = None

        def next(self):
            """
            (QueueAsLinkedList.Iterator) -> Object
            Returns the next element.
            """
            if self._position is None:
                self._position = self._container._list.head
            else:
                self._position = self._position.next
            if self._position is None:
                raise StopIteration
            return self._position.datum

    def __iter__(self):
        """
        (QueueAsLinkedList) -> QueueAsLinkedList.Iterator
        Returns an iterator for this queue.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (QueueAsLinkedList, QueueAsLinkedList) -> int

        Comparse this queue with the given queue.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "QueueAsLinkedList test program."
        print QueueAsLinkedList.main.__doc__
        queue2 = QueueAsLinkedList()
        Queue.test(queue2)
        return 0

if __name__ == "__main__":
    sys.exit(QueueAsLinkedList.main(*sys.argv))
