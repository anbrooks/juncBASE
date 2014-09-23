#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: naryTree.py,v $
#   $Revision: 1.23 $
#
#   $Id: naryTree.py,v 1.23 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the NaryTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.23 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object
from opus7.tree import Tree
from opus7.array import Array
from opus7.exception import *

#{
class NaryTree(Tree):
    """
    N-ary tree implemented using an array of._subtrees.
    """

#}@head


#{

    # ...
#}@tail

#{
    def __init__(self, *args):
        """
        (NaryTree, int [, Object]) -> None
        Constructs an N-ary tree.
        """
        super(NaryTree, self).__init__()
        if len(args) == 1:
            self._degree = args[0]
            self._key = None
            self._subtree = None
        elif len(args) == 2:
            self._degree = args[0]
            self._key = args[1]
            self._subtree = Array(self._degree)
            for i in xrange(self._degree):
                self._subtree[i] = NaryTree(self._degree)

    def purge(self):
        """
        (NaryTree) -> None
        Purges this N-ary tree.
        """
        self._key = None
        self._subtree = None
#}>a


    def getDegree(self):
        """
        (NaryTree) -> int
        Returns the degree of this N-ary tree.
        """
        if self.isEmpty:
            return 0
        else:
            return self._degree

    def getIsLeaf(self):
        """
        (NaryTree) -> bool
        Returns true if this N-ary tree is a leaf.
        """
        if self.isEmpty:
            return False
        for i in xrange(self._degree):
            if not self._subtree[i].isEmpty:
                return False
        return True

#{
    def getIsEmpty(self):
        """
        (NaryTree) -> bool
        Returns true if this N-ary tree is empty.
        """
        return self._key is None

    def getKey(self):
        """
        (NaryTree) -> Object
        Returns the._key of this N-ary tree node.
        """
        if self.isEmpty:
            raise StateError
        return self._key

    def attachKey(self, obj):
        """
        (NaryTree, Object) -> None
        Makes the given object the._key of this N-ary tree node.
        """
        if not self.isEmpty:
            raise StateError
        self._key = obj
        self._subtree = Array(self._degree)
        for i in xrange(self._degree):
            self._subtree[i] = NaryTree(self._degree)

    def detachKey(self):
        """
        (NaryTree) -> Object
        Detaches and returns the._key of this N-ary tree node.
        """
        if not isLeaf:
            raise StateError
        result = self._key
        self._key = None
        self._subtree = None
        return result
#}>b

#{
    def getSubtree(self, i):
        """
        (NaryTree, int) -> NaryTree
        Returns the specified._subtree of this N-ary tree.
        """
        if self.isEmpty:
            raise StateError
        return self._subtree[i]

    def attachSubtree(self, i, t):
        """
        (NaryTree, int, NaryTree) -> None
        Attaches the given tree as the specified._subtree of this N-ary tree.
        """
        if self.isEmpty or not self._subtree[i].isEmpty:
            raise StateError
        self._subtree[i] = t

    def detachSubtree(self, i):
        """
        (NaryTree, int) -> NaryTree
        Detaches and returns the specified._subtree of this N-ary tree.
        """
        if self.isEmpty:
            raise StateError
        result = self._subtree[i]
        self._subtree[i] = NaryTree(degree)
        return result
#}>c

    def _compareTo(self, obj):
        """
        (NaryTree, NaryTree) -> int

        Compares this N-ary tree with the given N-ary tree.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "NaryTree test program."
        print NaryTree.main.__doc__
        nt = NaryTree(3, 1)
        nt.attachSubtree(0, NaryTree(3, 2))
        nt.attachSubtree(1, NaryTree(3, 3))
        nt.attachSubtree(2, NaryTree(3, 4))
        Tree.test(nt)
        return 0

if __name__ == "__main__":
    sys.exit(NaryTree.main(*sys.argv))
