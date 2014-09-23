#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: mWayTree.py,v $
#   $Revision: 1.35 $
#
#   $Id: mWayTree.py,v 1.35 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the MWayTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.35 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.tree import Tree
from opus7.searchTree import SearchTree
from opus7.array import Array
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.stackAsLinkedList import StackAsLinkedList
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.prePostVisitor import PrePostVisitor
from opus7.exception import *

#{
class MWayTree(SearchTree):
    """
    M-way tree implemented using arrays.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, m):
        """
        (MWayTree, int) -> None
        Constructs an empty M-way tree.
        """
	assert m > 2
        super(MWayTree, self).__init__()
        self._key = Array(m - 1, 1)
        self._subtree = Array(m)

    def getM(self):
        """
        (MWayTree) -> int
        Returns the value of M for this M-way tree.
        """
        return len(self._subtree)

    m = property(
        fget = lambda self: self.getM())
#}>a

    def purge(self):
        """
        (MWayTree) -> None
        Purges this M-way tree.
        """
        for i in xrange(1, self._count + 1):
            self._key[i] = None
        for i in xrange(0, self._count + 1):
            self._subtree[i] = None
        self._count = 0

    def breadthFirstTraversal(self, visitor):
        """
        (MWayTree, Visitor) -> None
        Makes the given visitor do a breadth-first traversal
        of this M-way tree.
        """
        assert isinstance(visitor, Visitor)
        queue = QueueAsLinkedList()
        if not self.isEmpty:
            queue.enqueue(self)
        while not queue.isEmpty:
            head = queue.dequeue()
            for i in xrange(1, head.degree):
                visitor.visit(head.getKey(i))
            for i in xrange(0, head.degree):
                child = head.getSubtree(i)
                if not child.isEmpty:
                    queue.enqueue(child)

    def breadthFirstGenerator(self):
        """
        (BinaryTree) -> generator
        Yields the keys in this tree in breadth-first traversal order.
        """
        queue = QueueAsLinkedList()
        if not self.isEmpty:
            queue.enqueue(self)
        while not queue.isEmpty:
            head = queue.dequeue()
            for i in xrange(1, head.degree):
                yield head.getKey(i)
            for i in xrange(0, head.degree):
                child = head.getSubtree(i)
                if not child.isEmpty:
                    queue.enqueue(child)


    def getIsEmpty(self):
        """
        (MWayTree) -> bool
        Returns true if this M-way tree is empty.
        """
        return self._count == 0

    def getIsFull(self):
        """
        (MWayTree) -> bool
        Returns true if this M-way tree is full.
        """
        return self._count == self.m - 1

    def getIsLeaf(self):
        """
        (MWayTree) -> bool
        Returns true if this M-way tree is a leaf.
        """
        if self.isEmpty:
            return False
        for i in xrange(0, self._count + 1):
            if not self._subtree[i].isEmpty:
                return False
        return True

    def getDegree(self):
        """
        (MWayTree) -> int
        Returns the degree of this M-way tree node.
        """
        if self._count == 0:
            return 0
        else:
            return self._count + 1

    def getKey(self, *args):
        """
        (MWayTree, ...) -> Object
        Returns the specified key of this M-way tree node.
        """
        if self.isEmpty:
            raise StateError
	if len(args) == 0:
	    return self.getKey(0)
	elif len(args) == 1:
	    return self._key[args[0]]
	else:
	    raise ValueError

    def getSubtree(self, i):
        """
        (MwayTree, int) -> MWayTree
        Returns the specified subtree of this M-way tree node.
        """
        if self.isEmpty:
            raise StateError
        return self._subtree[i]

#{
    def depthFirstTraversal(self, visitor):
        """
        (MWayTree, PrePostVisitor) -> None
        Makes the given visitor do a depth-first traversal of this M-way tree.
        """
        assert isinstance(visitor, PrePostVisitor)
        if not self.isEmpty:
            for i in xrange(self._count + 2):
                if i > 1:
                    visitor.postVisit(self._key[i - 1])
                if i >= 1 and i <= self._count:
                    visitor.inVisit(self._key[i])
                if i < self._count:
                    visitor.preVisit(self._key[i + 1])
                if i <= self._count:
                    self._subtree[i].depthFirstTraversal(visitor)
#}>b

    def depthFirstGenerator(self, mode):
        """
        (BinaryTree) -> generator
        Yields the keys in this tree in depth-first traversal order.
        """
        if not self.isEmpty:
            for i in xrange(self._count + 2):
                if mode == self.POSTORDER and i > 1:
                    yield self._key[i - 1]
                if mode == self.INORDER and i >= 1 and i <= self._count:
                    yield self._key[i]
                if mode == self.PREORDER and i < self._count:
                    yield self._key[i + 1]
                if i <= self._count:
                    for obj in self._subtree[i].depthFirstGenerator(mode):
                        yield obj

    def getCount(self):
        """
        (MWayTree) -> int
        Returns the number of keys in this M-way tree node.
        """
        if self.isEmpty:
            return 0
        result = self._count
        for i in xrange(self._count + 1):
            result = result + self._subtree[i]._count
        return result

#{
    def find(self, obj):
        """
        (MWayTree, Object) -> Object
        Returns the object in this M-way tree that matches the given object.
        """
        if self.isEmpty:
            return None
        i = self._count
        while i > 0:
            diff = cmp(obj, self._key[i])
            if diff == 0:
                return self._key[i]
            if diff > 0:
                break
            i = i - 1
        return self._subtree[i].find(obj)
#}>c

    def __contains__(self, obj):
        """
        (MWayTree, Object) -> bool
        Returns true if the given object is in this M-way tree.
        """
        if self.isEmpty:
            return False
        i = self._count
        while i > 0:
            if self._key[i] is obj:
                return True
            if obj > self._key[i]:
                break
            i = i - 1
        return obj in self._subtree[i]

    def getMin(self):
        """
        (MWayTree) -> Object
        Returns the object in this M-way tree with the smallest value.
        """
        if self.isEmpty:
            return None
        elif self._subtree[0].isEmpty:
            return self._key[1]
        else:
            return self._subtree[0].min

    def getMax(self):
        """
        (MWayTree) -> Object
        Returns the object in this M-way tree with the largest value.
        """
        if self.isEmpty:
            return None
        elif self._subtree[self._count].isEmpty:
            return self._key[self._count]
        else:
            return self._subtree[self._count].max

#{
    def findIndex(self, obj):
        """
        (MWayTree, Object) -> int
        Returns the position of the specified object
        in the array of keys contained in this M-way tree node.
        Uses a binary search.
        """
        if self.isEmpty or obj < self._key[1]:
            return 0
        left = 1
        right = self._count
        while left < right:
            middle = (left + right + 1) / 2
            if obj < self._key[middle]:
                right = middle - 1
            else:
                left = middle
        return left

    def find(self, obj):
        """
        (MWayTree, Object):
        Returns the object in this M-way tree that matches the given object.
        """
        if self.isEmpty:
            return None
        index = self.findIndex(obj)
        if index != 0 and self._key[index] == obj:
            return self._key[index]
        else:
            return self._subtree[index].find(obj)
#}>d

#{
    def insert(self, obj):
        """
        (MWayTree, Object) -> None
        Inserts the given object into this M-way tree.
        """
        if self.isEmpty:
            self._subtree[0] = MWayTree(self.m)
            self._key[1] = obj
            self._subtree[1] = MWayTree(self.m)
            self._count = 1
        else:
            index = self.findIndex(obj)
            if index != 0 and self._key[index] == obj:
                raise ValueError
            if not self.isFull:
                i = self._count
                while i > index:
                    self._key[i + 1] = self._key[i]
                    self._subtree[i + 1] = self._subtree[i]
                    i = i - 1
                self._key[index + 1] = obj
                self._subtree[index + 1] = MWayTree(self.m)
                self._count = self._count + 1
            else:
                self._subtree[index].insert(obj)
#}>e

#{
    def withdraw(self, obj):
        """
        (MWayTree, Object) -> None
        Withdraws the given object from this M-way tree.
        """
        if self.isEmpty:
            raise KeyError
        index = self.findIndex(obj)
        if index != 0 and self._key[index] == obj:
            if not self._subtree[index - 1].isEmpty:
                max = self._subtree[index - 1].max
                self._key[index] = max
                self._subtree[index - 1].withdraw(max)
            elif not self._subtree[index].isEmpty:
                min = self._subtree[index].min
                self._key[index] = min
                self._subtree[index].withdraw(min)
            else:
                self._count = self._count - 1
                i = index
                while i <= self._count:
                    self._key[i] = self._key[i + 1]
                    self._subtree[i] = self._subtree[i + 1]
                    i = i + 1
                self._key[i] = None
                self._subtree[i] = None
                if self._count == 0:
                    self._subtree[0] = None
        else:
            self._subtree[index].withdraw(obj)
#}>f

    class Iterator(Iterator):
        """
        Enumerates the objects in an M-way tree.
        """

        def __init__(self, tree):
            """
            (MWayTree.Iterator, MWayTree) -> None
            Constructs an iterator for the given M-way tree.
            """
            super(MWayTree.Iterator, self).__init__(tree)
            self._position = 0
            self._stack = StackAsLinkedList()
        
        def next(self):
            """
            (MWayTree.Iterator) -> Object
            Returns the next object in this M-way tree.
            """
            if self._stack.isEmpty:
                if not self._container.isEmpty:
                    self._stack.push(self._container)
                self._position = 1
            else:
                self._position = self._position + 1
                top = self._stack.top
                if self._position == top.degree:
                    pop = self._stack.pop()
                    i = pop.degree - 1
                    while i >= 0:
                        subtree = pop.getSubtree(i)
                        if not subtree.isEmpty:
                            self._stack.push(subtree)
                        i = i - 1
                    self._position = 1
                if self._stack.isEmpty:
                    raise StopIteration
            return self._stack.top.getKey(self._position)

    def __iter__(self):
        """
        (MWayTree) -> MWayTree.Iterator
        Returns an iterator for this M-way tree.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (MWayTree, MWayTree) -> int

        Compares this M-way with the given M-way.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "MWayTree test program."
        print MWayTree.main.__doc__
        tree = MWayTree(3)
        SearchTree.test(tree)
        return 0

if __name__ == "__main__":
    sys.exit(MWayTree.main(*sys.argv))
