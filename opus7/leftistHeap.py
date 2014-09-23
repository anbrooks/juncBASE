#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: leftistHeap.py,v $
#   $Revision: 1.12 $
#
#   $Id: leftistHeap.py,v 1.12 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the LeftistHeap class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.binaryTree import BinaryTree
from opus7.priorityQueue import PriorityQueue
from opus7.mergeablePriorityQueue import MergeablePriorityQueue
from opus7.exception import *

#{
class LeftistHeap(BinaryTree, MergeablePriorityQueue):
    """
    Mergeable priority queue implemented as a leftist heap binary tree.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, *args):
        """
        (LeftistHeap [, key]) -> None
        Constructor.
        """
        if len(args) == 0:
            super(LeftistHeap, self).__init__()
            self._nullPathLength = 0
        else:
            super(LeftistHeap, self).__init__(
                args[0], LeftistHeap(), LeftistHeap())
            self._nullPathLength = 1
#}>a

    def swapContentsWith(self, heap):
        """
        (LeftistHeap, LeftistHeap) -> None
        Swaps the contents of this leftist heap with the given leftist heap.
        """
        tmp = self._key
        self._key = heap._key
        heap._key = tmp
        tmp = self._left
        self._left = heap._left
        heap._left = tmp
        tmp = self._right
        self._right = heap._right
        heap._right = tmp
        tmp = self._nullPathLength
        self._nullPathLength = heap._nullPathLength
        heap._nullPathLength = tmp

    def swapSubtrees(self):
        """
        (LeftistHeap) -> None
        Swaps the subtrees of this leftist heap.
        """
        tmp = self._left
        self._left = self._right
        self._right = tmp

#{
    def merge(self, queue):
        """
        (LeftistHeap, LeftistHeap) -> None
        Merges the contents of the given leftist heap with this leftist heap.
        """
        if self.isEmpty:
            self.swapContentsWith(queue)
        elif not queue.isEmpty:
            if self._key > queue._key:
                self.swapContentsWith(queue)
            self._right.merge(queue)
            if self._left._nullPathLength \
                    < self._right._nullPathLength:
                self.swapSubtrees()
            self._nullPathLength = 1 + min(
                self._left._nullPathLength,
                self._right._nullPathLength)
#}>b

#{
    def enqueue(self, obj):
        """
        (LeftistHeap, Object) -> None
        Enqueues the given object in this leftist heap.
        """
        self.merge(LeftistHeap(obj))
#}>c

#{
    def getMin(self):
        """
        (LeftistHeap) -> Object
        Returns the object in this leftist heap with the smallest value.
        """
        if self.isEmpty:
            raise ContainerEmpty
        return self._key
#}>d

#{
    def dequeueMin(self):
        """
        (LeftistHeap) -> Object
        Dequeues and returns the object in this leftist heap
        with the smallest value.
        """
        if self.isEmpty:
            raise ContainerEmpty
        result = self._key
        oldLeft = self._left
        oldRight = self._right
        self.purge()
        self.swapContentsWith(oldLeft)
        self.merge(oldRight)
        return result
#}>e

    @staticmethod
    def main(*argv):
        "LeftistHeap test program."
        print LeftistHeap.main.__doc__
        pqueue = LeftistHeap()
        PriorityQueue.test(pqueue)
        return 0

if __name__ == "__main__":
    sys.exit(LeftistHeap.main(*sys.argv))
