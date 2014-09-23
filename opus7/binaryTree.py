#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: binaryTree.py,v $
#   $Revision: 1.28 $
#
#   $Id: binaryTree.py,v 1.28 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BinaryTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.28 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.tree import Tree
from opus7.prePostVisitor import PrePostVisitor
from opus7.exception import *

#{
class BinaryTree(Tree):
    """
    Binary tree class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, *args):
        """
        (BinaryTree [, Object [, BinaryTree, BinaryTree]]) -> None
        Constructs a binary tree.
        """
        super(BinaryTree, self).__init__()
        if len(args) == 0:
            self._key = None
            self._left = None
            self._right = None
        elif len(args) == 1:
            self._key = args[0]
            self._left = BinaryTree()
            self._right = BinaryTree()
        elif len(args) == 3:
            self._key = args[0]
            self._left = args[1]
            self._right = args[2]
        else:
            raise ValueError

    def purge(self):
        """
        (BinaryTree) -> None
        Purges this binary tree.
        """
        self._key = None
        self._left = None
        self._right = None
#}>a

#{
    def getLeft(self):
        """
        (BinaryTree) -> BinaryTree
        Returns the left subtree of this binary tree.
        """
        if self.isEmpty:
            raise StateError
        return self._left

    left = property(
        fget = lambda self: self.getLeft())

    def getRight(self):
        """
        (BinaryTree) -> BinaryTree
        Returns the right subtree of this binary tree.
        """
        if self.isEmpty:
            raise StateError
        return self._right

    right = property(
        fget = lambda self: self.getRight())
#}>b

    def getIsLeaf(self):
        """
        (BinaryTree) -> bool
        Returns true if this binary tree node is a leaf.
        """
        return not self.isEmpty and self._left.isEmpty and \
            self._right.isEmpty

    def getDegree(self):
        """
        (BinaryTree) -> int
        Returns the degree of this binary tree node.
        """
        if self.isEmpty:
            return 0
        else:
            return 2

    def getIsEmpty(self):
        """
        (BinaryTree) -> bool
        Returns true if this binary tree is empty.
        """
        return self._key is None

    def getKey(self):
        """
        (BinaryTree) -> Object
        Returns the key in this binary tree node.
        """
        if self.isEmpty:
            raise ValueError
        return self._key

    def getSubtree(self, i):
        """
        (BinaryTree, int) -> Object
        Returns the specified subtree of this binary tree.
        """
        if i == 0:
            return self._left
        elif i == 1:
            return self._right
        else:
            raise IndexError

    def attachKey(self, obj):
        """
        (BinaryTree, Object) -> None
        Makes the given object the key of this binary tree node.
        """
        if not self.isEmpty:
            raise StateError
        self._key = obj
        self._left = BinaryTree()
        self._right = BinaryTree()

    def detachKey(self):
        """
        (BinaryTree) -> None
        Detaches and returns the key in this binary tree node.
        """
        if not self.isLeaf:
            raise StateError
        result = self._key
        self._key = None
        self._left = None
        self._right = None
        return result

    def attachLeft(self, t):
        """
        (BinaryTree, Binary Tree) -> None
        Attaches the given binary tree as the left subtree of this binary tree.
        """
        if self.isEmpty or not self._left.isEmpty:
            raise StateError
        self._left = t

    def detachLeft(self):
        """
        (BinaryTree) -> BinaryTree
        Detaches and returns the left subtree of this binary tree.
        """
        if self.isEmpty:
            raise StateError
        result = self._left
        self._left = BinaryTree()
        return result

    def attachRight(self, t):
        """
        (BinaryTree, Binary Tree) -> None
        Attaches the given binary tree as the right subtree of this binary tree.
        """
        if self.isEmpty or not self._right.isEmpty:
            raise StateError
        self._right = t

    def detachRight(self):
        """
        (BinaryTree) -> BinaryTree
        Detaches and returns the right subtree of this binary tree.
        """
        if self.isEmpty:
            raise StateErorr
        result = self._right
        self._right = BinaryTree()
        return result

#{
    def depthFirstTraversal(self, visitor):
        """
        (BinaryTree, PrePostVisitor) -> None
        Makes the given visitor do a depth-first traversal of this binary tree.
        """
        assert isinstance(visitor, PrePostVisitor)
        if not self.isEmpty:
            visitor.preVisit(self.key)
            self.left.depthFirstTraversal(visitor)
            visitor.inVisit(self.key)
            self.right.depthFirstTraversal(visitor)
            visitor.postVisit(self.key)
#}>c

    def depthFirstGenerator(self, mode):
        """
        (BinaryTree) -> generator
        Yields the keys in this tree in depth-first traversal order.
        """
        if not self.isEmpty:
            if mode == self.PREORDER:
                yield self.key
            for obj in self.left.depthFirstGenerator(mode):
                yield obj
            if mode == self.INORDER:
                yield self.key
            for obj in self.right.depthFirstGenerator(mode):
                yield obj
            if mode == self.POSTORDER:
                yield self.key

#{
    def _compareTo(self, bt):
        """
        (BinaryTree, BinaryTree) -> int
        Compares this binary tree to the given binary tree.
        """
        assert isinstance(self, bt.__class__)
        if self.isEmpty:
            if bt.isEmpty:
                return 0
            else:
                return -1
        elif bt.isEmpty:
            return 1
        else:
            result = cmp(self._key, bt._key)
            if result == 0:
                result = cmp(self._left, bt._left)
            if result == 0:
                result = cmp(self._right, bt._right)
            return result
#}>d

    @staticmethod
    def main(*argv):
        "BinaryTree test program."
        print BinaryTree.main.__doc__
        bt = BinaryTree(4)
        bt.attachLeft(BinaryTree(2))
        bt.attachRight(BinaryTree(6))
        Tree.test(bt)
        return 0

if __name__ == "__main__":
    sys.exit(BinaryTree.main(*sys.argv))
