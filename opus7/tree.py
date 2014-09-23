#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: tree.py,v $
#   $Revision: 1.34 $
#
#   $Id: tree.py,v 1.34 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the Tree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.34 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.container import Container
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.stackAsLinkedList import StackAsLinkedList
from opus7.prePostVisitor import PrePostVisitor
from opus7.preOrder import PreOrder
from opus7.inOrder import InOrder
from opus7.postOrder import PostOrder
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.printingVisitor import PrintingVisitor

#{
class Tree(Container):
    """
    Base class from which all tree classes are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    @abstractmethod
    def getKey(self):
        pass

    key = property(
        fget = lambda self: self.getKey())

    @abstractmethod
    def getSubtree(self, i):
        pass

    @abstractmethod
    def getIsLeaf(self):
        pass

    isLeaf = property(
        fget = lambda self: self.getIsLeaf())

    @abstractmethod
    def getDegree(self):
        pass

    degree = property(
        fget = lambda self: self.getDegree())

    @abstractmethod
    def getHeight(self):
        pass

    height = property(
        fget = lambda self: self.getHeight())
#}>a

    def __init__(self):
        """
        (Tree) -> None
        Constructor.
        """
        super(Tree, self).__init__()

#{
    def depthFirstTraversal(self, visitor):
        """
        (Tree, PrePostVisitor) -> None
        Makes the given visitor do a depth-first traversal of this tree.
        """
        assert isinstance(visitor, PrePostVisitor)
        if not self.isEmpty and not visitor.isDone:
            visitor.preVisit(self.key)
            for i in xrange(self.degree):
                self.getSubtree(i).depthFirstTraversal(visitor)
            visitor.postVisit(self.key)
#}>b

    PREORDER = -1
    INORDER = 0
    POSTORDER = +1

    def depthFirstGenerator(self, mode):
        """
        (Tree) -> generator
        Yields the keys in this tree in depth-first traversal order.
        """
        if not self.isEmpty:
            if mode == self.PREORDER:
                yield self.key
            for i in xrange(self.degree):
                for obj in self.getSubtree(i).depthFirstGenerator(mode):
                    yield obj
            if mode == self.POSTORDER:
                yield self.key

#{
    def breadthFirstTraversal(self, visitor):
        """
        (Tree, Visitor) -> None
        Makes the given visitor do a breadth-first traversal of this tree.
        """
        assert isinstance(visitor, Visitor)
        queue = QueueAsLinkedList()
        if not self.isEmpty:
            queue.enqueue(self)
        while not queue.isEmpty and not visitor.isDone:
            head = queue.dequeue()
            visitor.visit(head.key)
            for i in xrange(head.degree):
                child = head.getSubtree(i)
                if not child.isEmpty:
                    queue.enqueue(child)
#}>c

    def breadthFirstGenerator(self):
        """
        (Tree) -> generator
        Yields the keys in this tree in depth-first traversal order.
        """
        queue = QueueAsLinkedList()
        if not self.isEmpty:
            queue.enqueue(self)
        while not queue.isEmpty:
            head = queue.dequeue()
            yield head.key
            for i in xrange(head.degree):
                child = head.getSubtree(i)
                if not child.isEmpty:
                    queue.enqueue(child)

#{
    def accept(self, visitor):
        """
        (Tree) -> Visitor
        Makes the given visitor visit the nodes of this tree.
        """
        assert isinstance(visitor, Visitor)
        self.depthFirstTraversal(PreOrder(visitor))
#}>d

    def getHeight(self):
        """
        (Tree) -> int
        Returns the height of this tree.
        """
        if self.isEmpty:
            return -1
        height = -1
        for i in xrange(self.degree):
            height = max(height, self.getSubtree(i).height)
        return height + 1

    def getCount(self):
        """
        (Tree) -> int
        Returns the number of nodes in this tree.
        """
        if self.IsEmpty():
            return 0
        result = 1
        for i in xrange(self.degree):
            result = result + self.getSubtree(i).count
        return result

#{
    class Iterator(Iterator):
        """
        Enumerates the nodes of a tree.
        """

        def __init__(self, tree):
            """
            (Tree.Iterator, Tree) -> None
            Constructs an iterator for the given tree.
            """
            super(Tree.Iterator, self).__init__(tree)
            self._stack = StackAsLinkedList()
	    if not tree.isEmpty:
		self._stack.push(tree)

        def next(self):
            """
            (Tree.Iterator) -> Object
            Returns the next node in the tree.
            """
            if self._stack.isEmpty:
		raise StopIteration
	    top = self._stack.pop()
	    i = top.degree - 1
	    while i >= 0:
		subtree = top.getSubtree(i)
		if not subtree.isEmpty:
		    self._stack.push(subtree)
		i -= 1
            return top.key

    def __iter__(self):
        """
        (Tree) -> Tree.Iterator
        Returns an interator for this tree.
        """
        return Tree.Iterator(self)
#}>e

    @staticmethod
    def test(tree):
        "Tree test program."
        print Tree.test.__doc__
        visitor = PrintingVisitor()
        print tree
        print "Breadth-First traversal"
        tree.breadthFirstTraversal(visitor)
        visitor.finish()
        print "Preorder traversal"
        tree.depthFirstTraversal(PreOrder(visitor))
        visitor.finish()
        print "Inorder traversal"
        tree.depthFirstTraversal(InOrder(visitor))
        visitor.finish()
        print "Postorder traversal"
        tree.depthFirstTraversal(PostOrder(visitor))
        visitor.finish()
        print "Using iterator"
        for i in tree:
            print i
        print "Using depth-first generator (preorder)"
        for i in tree.depthFirstGenerator(Tree.PREORDER):
            print i
        print "Using breadth-first generator"
        for i in tree.breadthFirstGenerator():
            print i
