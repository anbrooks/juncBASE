#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: dequeAsLinkedList.py,v $
#   $Revision: 1.21 $
#
#   $Id: dequeAsLinkedList.py,v 1.21 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the DequeAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.21 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.queue import Queue
from opus7.deque import Deque
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.deque import Deque
from opus7.exception import *

#{
class DequeAsLinkedList(QueueAsLinkedList, Deque):
    """
    Deque implemented using a linked list.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (DequeAsLinkedList) -> None
        Constructs a deque.
        """
        super(DequeAsLinkedList, self).__init__()

    def enqueueHead(self, obj):
        """
        (DequeAsLinkedList, Object) -> None
        Enqueues the given object at the head of this deque.
        """
        self._list.prepend(obj)
        self._count += 1
#}>a

#{
    def getTail(self):
        """
        (DequeAsLinkedList) -> Object
        Returns the object at the tail of this deque.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._list.last

    def dequeueTail(self):
        """
        (DequeAsLinkedList) -> Object
        Dequeues the object at the tail of this deque.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._list.last
        self._list.extract(result)
        self._count -= 1
        return result
#}>b

    @staticmethod
    def main(*argv):
        "DequeAsLinkedList test program."
        print DequeAsLinkedList.main.__doc__
        deque2 = DequeAsLinkedList()
        Queue.test(deque2)
        Deque.test(deque2)
        return 0

if __name__ == "__main__":
    sys.exit(DequeAsLinkedList.main(*sys.argv))
