#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: bTree.py,v $
#   $Revision: 1.21 $
#
#   $Id: bTree.py,v 1.21 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.21 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.searchTree import SearchTree
from opus7.mWayTree import MWayTree

#{
class BTree(MWayTree):
    """
    B-tree class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, m):
        """
        (BTree, int) -> None
        Constructs an empty B-tree with the given value of M.
        """
        super(BTree, self).__init__(m)
        self._parent = None

    def attachSubtree(self, i, t):
        """
        (BTree, int, Btree):
        Attaches the given B-tree as the specified subtree of this B-tree.
        """
        self._subtree[i] = t
        t._parent = self
#}>a

#{
    def insert(self, obj):
        """
        (Btree, Object) -> None
        Inserts the given object into this B-tree.
        """
        if self.isEmpty:
            if self._parent is None:
                self.attachSubtree(0, BTree(self.m))
                self._key[1] = obj
                self.attachSubtree(1, BTree(self.m))
                self._count = 1
            else:
                self._parent.insertUp(obj, BTree(self.m))
        else:
            index = self.findIndex(obj)
            if index != 0 and self._key == obj:
                raise KeyError
            self._subtree[index].insert(obj)
#}>b

#{
    def insertUp(self, obj, child):
        """
        (BTree, Object, Btree) -> None
        Inserts the given (Object, Btree) pair into this B-tree.
        """
        index = self.findIndex(obj)
        if not self.isFull:
	    self.insertPair(index + 1, obj, child)
            self._count = self._count + 1
        else:
            (extraKey, extraTree) = self.insertPair(
		index + 1, obj, child)
            if self._parent is None:
                left = BTree(self.m)
                right = BTree(self.m)
                left.attachLeftHalfOf(self)
                right.attachRightHalfOf(self)
                right.insertUp(extraKey, extraTree)
                self.attachSubtree(0, left)
                self._key[1] = self._key[(self.m + 1)/2]
                self.attachSubtree(1, right)
                self._count = 1
            else:
                self._count = (self.m + 1)/2 - 1
                right = BTree(self.m)
                right.attachRightHalfOf(self)
                right.insertUp(extraKey, extraTree)
                self._parent.insertUp(
                    self._key[(self.m + 1)/2], right)
#}>c

    def insertPair(self, index, obj, child):
        """
        (BTree, int, Object, BTree) -> (Object, BTree)
        Inserts the given Object at the specified index
        in the key array of this B-tree node
        and returns any leftover object.
        """
        if index == self.m:
            return (obj, child)
        result = (self._key[self.m - 1], self._subtree[self.m - 1])
        i = self.m - 1
        while i > index:
            self._key[i] = self._key[i - 1]
            self._subtree[i] = self._subtree[i - 1]
            i = i - 1
        self._key[index] = obj
        self._subtree[index] = child
        child._parent = self
        return result

    def attachLeftHalfOf(self, btree):
        """
        (BTree, Btree) -> None
        Attaches the left half of this B-tree to the given B-tree.
        """
        self._count = (self.m + 1)/2 - 1
        self.attachSubtree(0, btree._subtree[0])
        for i in xrange(1, self._count + 1):
            self._key[i] = btree._key[i]
            self.attachSubtree(i, btree._subtree[i])

    def attachRightHalfOf(self, btree):
        """
        (BTree, Btree) -> None
        Attaches the right half of this B-tree to the given B-tree.
        """
        self._count = self.m - (self.m + 1)/2 - 1
        j = (self.m + 1)/2
        self.attachSubtree(0, btree._subtree[j])
        j = j + 1
        for i in xrange(1, self._count + 1):
            self._key[i] = btree._key[j]
            self.attachSubtree(i, btree._subtree[j])

    def withdraw(self, obj):
        """
        (Btree, Object) -> Object

        Withdraws the given object from this B-tree.
        """
        raise NotImplementedError
    
    @staticmethod
    def main(*argv):
        "BTree test program."
        print BTree.main.__doc__
        tree = BTree(3)
        SearchTree.test(tree)
        return 0

if __name__ == "__main__":
    sys.exit(BTree.main(*sys.argv))
