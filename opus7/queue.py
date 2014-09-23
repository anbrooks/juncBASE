#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: queue.py,v $
#   $Revision: 1.27 $
#
#   $Id: queue.py,v 1.27 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Queue class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.27 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.container import Container

#{
class Queue(Container):
    """
    Base class from which all queue classes are derived.
    """

    def __init__(self):
        """
        (Queue) -> None
        Constructor.
        """
        super(Queue, self).__init__()

    @abstractmethod
    def getHead(self): pass

    head = property(
        fget = lambda self: self.getHead())

    @abstractmethod
    def enqueue(self, obj): pass

    @abstractmethod
    def dequeue(self): pass
#}>a

    @staticmethod
    def test(queue):
        "Queue test program."
        print Queue.test.__doc__
        for i in xrange(6):
            if queue.isFull:
                break
            queue.enqueue(i)
        print queue
        print queue.head
        while not queue.isEmpty:
            obj = queue.dequeue()
            print obj
