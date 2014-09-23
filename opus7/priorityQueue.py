#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: priorityQueue.py,v $
#   $Revision: 1.16 $
#
#   $Id: priorityQueue.py,v 1.16 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PriorityQueue class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.container import Container

#{
class PriorityQueue(Container):
    """
    Base class from which all priority queue classes are derived.
    """

    def __init__(self):
        """
        (PriorityQueue) -> None
        Constructor.
        """
        super(PriorityQueue, self).__init__()

    @abstractmethod
    def enqueue(self, obj): pass

    @abstractmethod
    def getMin(self): pass

    min = property(
        fget = lambda self: self.getMin())

    @abstractmethod
    def dequeueMin(self): pass
#}>a

    @staticmethod
    def test(pqueue):
        "PriorityQueue test program."
        print PriorityQueue.test.__doc__
        print pqueue
        pqueue.enqueue(3)
        pqueue.enqueue(1)
        pqueue.enqueue(4)
        pqueue.enqueue(1)
        pqueue.enqueue(5)
        pqueue.enqueue(9)
        pqueue.enqueue(2)
        pqueue.enqueue(6)
        pqueue.enqueue(5)
        pqueue.enqueue(4)
        print pqueue
        while not pqueue.isEmpty:
            obj = pqueue.dequeueMin()
            print obj

