#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: dequeAsArray.py,v $
#   $Revision: 1.23 $
#
#   $Id: dequeAsArray.py,v 1.23 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the DequeAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.23 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.queue import Queue
from opus7.deque import Deque
from opus7.queueAsArray import QueueAsArray
from opus7.exception import *

#{
class DequeAsArray(QueueAsArray, Deque):
    """
    Deque implemented using an array.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size = 0):
        """
        (DequeAsArray, int) -> None
        Constructs a deque of the given size.
        """
        super(DequeAsArray, self).__init__(size)

    def enqueueHead(self, obj):
        """
        (DequeAsArray, Object) -> None
        Enqueues the given object at the head of this deque.
        """
        if self._count == len(self._array):
            raise ContainerFull
        if self._head == 0:
            self._head = len(self._array) - 1
        else:
            self._head = self._head - 1
        self._array[self._head] = obj
        self._count += 1
#}>a

#{
    def getTail(self):
        """
        (DequeAsArray) -> Object
        Returns the object at the tail of this deque.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._array[self._tail]

    def dequeueTail(self):
        """
        (DequeAsArray) -> Object
        Dequeues the objectc at the tail of this deque.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._array[self._tail]
        self._array[self._tail] = None
        if self._tail == 0:
            self._tail = len(self._array) - 1
        else:
            self._tail = self._tail - 1
        self._count -= 1
        return result
#}>b

    @staticmethod
    def main(*argv):
        "DequeAsArray test program."
        print DequeAsArray.main.__doc__
        deque1 = DequeAsArray(5)
        Queue.test(deque1)
        Deque.test(deque1)
        return 0

if __name__ == "__main__":
    sys.exit(DequeAsArray.main(*sys.argv))
