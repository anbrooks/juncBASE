#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: mergeablePriorityQueue.py,v $
#   $Revision: 1.11 $
#
#   $Id: mergeablePriorityQueue.py,v 1.11 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the MergeablePriorityQueue class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.priorityQueue import PriorityQueue

#{
class MergeablePriorityQueue(PriorityQueue):
    """
    A mergeable priority queue.
    """

    def __init__(self):
        """
        (MergeablePriorityQueue) -> None
        Constructor.
        """
        super(MergeablePriorityQueue, self).__init__()

    @abstractmethod
    def merge(self, queue): pass
#}>a

    def test(pqueue):
        "MergeablePriorityQueue test program."
        print MergeablePriorityQueue.test.__doc__
