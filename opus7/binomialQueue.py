#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: binomialQueue.py,v $
#   $Revision: 1.19 $
#
#   $Id: binomialQueue.py,v 1.19 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BinomialQueue class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.19 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.priorityQueue import PriorityQueue
from opus7.mergeablePriorityQueue import MergeablePriorityQueue
from opus7.linkedList import LinkedList
from opus7.generalTree import GeneralTree
from opus7.preOrder import PreOrder
from opus7.exception import *

#{
class BinomialQueue(MergeablePriorityQueue):
    """
    Mergeable priority queue implemented as a binomial queue.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, *args):
        """
        (BinomialQueue, ...) -> None
        Constructor.
        """
        super(BinomialQueue, self).__init__()
        self._treeList = LinkedList()
	if len(args) == 0:
	    pass
	elif len(args) == 1:
	    assert isinstance(args[0], self.BinomialTree)
	    self._treeList.append(args[0])
	else:
	    raise ValueError
#}>a


#{
    def addTree(self, tree):
        """
        (BinomialQueue, BinomialQueue.BinomialTree) -> None
        Adds the given binomial tree to this binomial queue.
        """
        self._treeList.append(tree)
        self._count += tree.count

    def removeTree(self, tree):
        """
        (BinomialQueue, BinomialQueue.BinomialTree) -> None
        Removes the given binomial tree from this binomial queue.
        """
        self._treeList.extract(tree)
        self._count -= tree.count
#}>b

    def purge(self):
        """
        (BinomialQueue) -> None
        Purges this binomial queue.
        """
        self._treeList = LinkedList()
        self._count = 0

    def accept(self, visitor):
        """
        (BinomialQueue, Visitor) -> None
        Makes the given visitor visit the elements of this binomial queue.
        """
        assert isinstance(visitor, Visitor)
        ptr = self._treeList.head
        while ptr is not None:
            tree = ptr.datum
            tree.depthFirstTraversal(PreOrder(visitor))
            ptr = ptr.next

#{
    def getMinTree(self):
        """
        (BinomialQueue) -> BinomialQueue.BinomialTree
        Returns the binomial tree in this binomial queue
        with the smallest root.
        """
        minTree = None
        ptr = self._treeList.head
        while ptr is not None:
            tree = ptr.datum
            if minTree is None or tree.key < minTree.key:
                minTree = tree
            ptr = ptr.next
        return minTree

    minTree = property(
        fget = lambda self: self.getMinTree())

    def getMin(self):
        """
        (BinomialQueue) -> Object
        Returns the object in this binomial queue with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self.minTree.key
#}>c

#{
    def merge(self, queue):
        """
        (BinomialQueue, BinomialQueue) -> None
        Merges the contents of the given binomial queue
        with this binomial queue.
        """
        oldList = self._treeList
        self._treeList = LinkedList()
        self._count = 0
        p = oldList.head
        q = queue._treeList.head
        carry = None
        i = 0
        while p is not None or q is not None \
                or carry is not None:
            a = None
            if p is not None:
                tree = p.datum
                if tree.degree == i:
                    a = tree
                    p = p.next
            b = None
            if q is not None:
                tree = q.datum
                if tree.degree == i:
                    b = tree
                    q = q.next
            (sum, carry) = BinomialQueue.fullAdder(a, b, carry)
            if sum is not None:
                self.addTree(sum)
            i += 1
        queue.purge()
#}>d


#{
    @staticmethod
    def fullAdder(a, b, c):
	"""
        (BinomialTree, BinomialTree, BinomialTree) ->
	    (BinomialTree, BinomialTree)
        Returns the (sum, carry) of the given binomial trees.
	"""
	if a is None:
	    if b is None:
		if c is None:
		    return (None, None)
		else:
		    return (c, None)
	    else:
		if c is None:
		    return (b, None)
		else:
		    return (None, b.add(c))
	else:
	    if b is None:
		if c is None:
		    return (a, None)
		else:
		    return (None, a.add(c))
	    else:
		if c is None:
		    return (None, a.add(b))
		else:
		    return (c, a.add(b))
#}>e

#{
    def enqueue(self, obj):
        """
        (BinomialQueue, Object) -> None
        Enqueues the given object in this binomial queue.
        """
        self.merge(BinomialQueue(
	    BinomialQueue.BinomialTree(obj)))
#}>f

#{
    def dequeueMin(self):
        """
        (BinomialQueue) -> Object
        Dequeues and returns the object in this binomial queue
        with the smallest value.
        """
        if self._count == 0:
            raise ContainerEmpty
        minTree = self.minTree
        self.removeTree(minTree)
        queue = BinomialQueue()
        while minTree.degree > 0:
            child = minTree.getSubtree(0)
            minTree.detachSubtree(child)
            queue.addTree(child)
        self.merge(queue)
        return minTree.key
#}>g

    def __str__(self):
        """
        (BinomialQueue) -> str
        Returns a string representation of this binomial queue.
        """
        result = self.__class__.__name__ + " {\n"
        ptr = self._treeList.head
        while ptr is not None:
            result = result + str(ptr.datum) + "\n"
            ptr = ptr.next
        result = result + "}"
        return result

#{
    class BinomialTree(GeneralTree):
        """
        A binomial tree implemented as a general tree.
        """

#}@head

#{

        # ...
#}@tail

#{
        def __init__(self, key):
            """
            (BinomialQueue.BinomialTree, Object) -> None
            Constructor.
            """
            super(BinomialQueue.BinomialTree, self).__init__(key)
#}>h

        def getCount(self):
            """
            (BinomialQueue.BinomialTree) -> int
            Returns the number of objects in this binomial tree.
            """
            return 1 << self._degree

        def swapContentsWith(self, tree):
            """
            (BinomialQueue.BinomialTree, BinomialQueue.BinomialTree) -> None
            Swaps the contents of this binomial tree
            with the given binomial tree.
            """
            tmp = self._key
            self._key = tree._key
            tree._key = tmp
            tmp = self._list
            self._list = tree._list
            tree._list = tmp
            tmp = self._degree
            self._degree = tree._degree
            tree._degree = tmp

#{
        def add(self, tree):
            """
            (BinomialQueue.BinomialTree, BinomialQueue.BinomialTree) ->
                BinomialQueue.BinomialTree
            Adds this binomail tree and the given binomial tree.
            """
            if self._degree != tree._degree:
                raise ValueError
            if self._key > tree._key:
                self.swapContentsWith(tree)
            self.attachSubtree(tree)
            return self
#}>i

    def _compareTo(self, obj):
        """
        (BinomialQueue, BinomialQueue) -> int

        Compares this binomial queue with the given binomial queue.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    def __iter__(self):
        """
        (BinomialQueue) -> iterator

        Returns an iterator that enumerates the elements of this binomial queue.
        """
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "BinomialQueue test program."
        print BinomialQueue.main.__doc__
        pqueue = BinomialQueue()
        PriorityQueue.test(pqueue)
        return 0

if __name__ == "__main__":
    sys.exit(BinomialQueue.main(*sys.argv))
