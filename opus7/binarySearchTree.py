#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: binarySearchTree.py,v $
#   $Revision: 1.20 $
#
#   $Id: binarySearchTree.py,v 1.20 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BinarySearchTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.20 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.binaryTree import BinaryTree
from opus7.searchTree import SearchTree
from opus7.exception import *

#{
class BinarySearchTree(BinaryTree, SearchTree):
    """
    Binary search tree class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (BinarySearchTree) -> None
        Constructs an empty binary search tree.
        """
        super(BinarySearchTree, self).__init__()
#}>a

    def __contains__(self, obj):
        """
        (BinarySearchTree, Object) -> bool
        Returns true if the given object is in this binary search tree.
        """
        if self.isEmpty:
            return False
        elif self._key is obj:
            return True
        elif obj < self._key:
            return obj in self._left
        elif obj > self._key:
            return obj in self._right
        else:
            return False

#{
    def find(self, obj):
        """
        (BinarySearchTree, Object) -> Object
        Returns the object in this binary search tree
        that matches the given object.
        """
        if self.isEmpty:
            return None
        diff = cmp(obj, self._key)
        if diff == 0:
            return self._key
        elif diff < 0:
            return self._left.find(obj)
        elif diff > 0:
            return self._right.find(obj)

    def getMin(self):
        """
        (BinarySearchTree) -> Object
        Returns the object in this binary search tree with the smallest value.
        """
        if self.isEmpty:
            return None
        elif self._left.isEmpty:
            return self._key
        else:
            return self._left.min
#}>b

    def getMax(self):
        """
        (BinarySearchTree) -> Object
        Returns the object in this binary search tree with the largest value.
        """
        if self.isEmpty:
            return None
        elif self._right.isEmpty:
            return self._key
        else:
            return self._right.max

#{
    def insert(self, obj):
        """
        (BinarySearchTree, Object) -> None
        Inserts the given object into this binary search tree.
        """
        if self.isEmpty:
            self.attachKey(obj)
        else:
            diff = cmp(obj, self._key)
            if diff == 0:
                raise ValueError
            elif diff < 0:
                self._left.insert(obj)
            elif diff > 0:
                self._right.insert(obj)
        self.balance()

    def attachKey(self, obj):
        """
        (BinarySearchTree, Object) -> None
        Attaches the given object to the root node of this binary search tree.
        """
        if not self.isEmpty:
            raise StateError
        self._key = obj
        self._left = BinarySearchTree()
        self._right = BinarySearchTree()

    def balance(self):
        """
        (BinarySearchTree) -> None
        Balances this binary search tree.
        """
        pass
#}>c

#{
    def withdraw(self, obj):
        """
        (BinarySearchTree, Object) -> None
        Withdraws the given object from this binary search tree.
        """
        if self.isEmpty:
            raise KeyError
        diff = cmp(obj, self._key)
        if diff == 0:
            if not self._left.isEmpty:
                max = self._left.max
                self._key = max
                self._left.withdraw(max)
            elif not self._right.isEmpty:
                min = self._right.min
                self._key = min
                self._right.withdraw(min)
            else:
                self.detachKey()
        elif diff < 0:
            self._left.withdraw(obj)
        elif diff > 0:
            self._right.withdraw(obj)
        self.balance()
#}>d

    @staticmethod
    def main(*argv):
        "BinarySearchTree test program."
        print BinarySearchTree.main.__doc__
        tree = BinarySearchTree()
        SearchTree.test(tree)
        return 0

if __name__ == "__main__":
    sys.exit(BinarySearchTree.main(*sys.argv))
