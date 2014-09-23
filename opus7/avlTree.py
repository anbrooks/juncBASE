#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: avlTree.py,v $
#   $Revision: 1.19 $
#
#   $Id: avlTree.py,v 1.19 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the AVLTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.19 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.searchTree import SearchTree
from opus7.binarySearchTree import BinarySearchTree
from opus7.exception import *

#{
class AVLTree(BinarySearchTree):
    """
    AVL tree class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (AVLTree) -> None
        Constructs an empty AVL tree.
        """
        super(AVLTree, self).__init__()
        self._height = -1
#}>a

#{
    def getHeight(self):
        """
        (AVLTree) -> int
        Returns the height of this AVL tree.
        """
        return self._height

    def adjustHeight(self):
        """
        (AVLTree) -> None
        Adjusts the height of this AVL tree.
        """
        if self.isEmpty:
            self._height = -1
        else:
            self._height = 1 + max(
                self._left._height, self._right._height)

    def getBalanceFactor(self):
        """
        (AVLTree) -> int
        Returns the balance factor of this AVL tree node.
        """
        if self.isEmpty:
            return 0
        else:
            return self._left._height - self._right._height

    balanceFactor = property(
        fget = lambda self: self.getBalanceFactor())
#}>b

#{
    def doLLRotation(self):
        """
        (AVLTree) -> None
        Does an LL rotation at this AVL tree node.
        """
        if self.isEmpty:
            raise StateError
        tmp = self._right
        self._right = self._left
        self._left = self._right._left
        self._right._left = self._right._right
        self._right._right = tmp

        tmp = self._key
        self._key = self._right._key
        self._right._key = tmp

        self._right.adjustHeight()
        self.adjustHeight()
#}>c

    def doRRRotation(self):
        """
        (AVLTree) -> None
        Does an RR rotation at this AVL tree node.
        """
        if self.isEmpty:
            raise StateError
        tmp = self._left
        self._left = self._right
        self._right = self._left._right
        self._left._right = self._left._left
        self._left._left = tmp

        tmp = self._key
        self._key = self._left._key
        self._left._key = tmp

        self._left.adjustHeight()
        self.adjustHeight()

#{
    def doLRRotation(self):
        """
        (AVLTree) -> None
        Does an LR rotation at this AVL tree node.
        """
        if self.isEmpty:
            raise StateError
        self._left.doRRRotation()
        self.doLLRotation()
#}>d

    def doRLRotation(self):
        """
        (AVLTree) -> None
        Does an RL rotation at this AVL tree node.
        """
        if self.isEmpty:
            raise StateError
        self._right.doLLRotation()
        self.doRRRotation()

#{
    def balance(self):
        """
        (AVLTree) -> None
        Balances this AVL tree node.
        """
        self.adjustHeight()
        if self.balanceFactor > 1:
            if self._left.balanceFactor > 0:
                self.doLLRotation()
            else:
                self.doLRRotation()
        elif self.balanceFactor < -1:
            if self._right.balanceFactor < 0:
                self.doRRRotation()
            else:
                self.doRLRotation()
#}>e

#{
    def attachKey(self, obj):
        """
        (AVLTree, Object) -> None
        Attaches the given object to this AVL tree node.
        """
        if not self.isEmpty:
            raise StateError
        self._key = obj
        self._left = AVLTree()
        self._right = AVLTree()
        self._height = 0
#}>f

#{
    def detachKey(self):
        """
        (AVLTree) -> Object
        Detaches and returns the key in this AVL tree node.
        """
        self._height = -1
        return super(AVLTree, self).detachKey()
#}>g
    
    @staticmethod
    def main(*argv):
        "AVLTree test program."
        print AVLTree.main.__doc__
        tree = AVLTree()
        SearchTree.test(tree)
        return 0

if __name__ == "__main__":
    sys.exit(AVLTree.main(*sys.argv))
