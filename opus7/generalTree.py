#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: generalTree.py,v $
#   $Revision: 1.24 $
#
#   $Id: generalTree.py,v 1.24 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the GeneralTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.24 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.tree import Tree
from opus7.linkedList import LinkedList

#{
class GeneralTree(Tree):
    """
    A general tree implemented using a linked list of subtrees.
    """

#}@head

#{
    # ...
#}@tail

#{
    def __init__(self, key):
        """
        (GeneralTree, Object) -> None
        Constructs a general tree with the given object at its root.
        """
        super(GeneralTree, self).__init__()
        self._key = key
        self._degree = 0
        self._list = LinkedList()

    def purge(self):
        """
        (GeneralTree) -> None
        Purges this general tree.
        """
        self._list.purge()
        self._degree = 0
#}>a

    def getIsEmpty(self):
        """
        (GeneralTree) -> bool
        Returns false always.
        """
        return False

    def getIsLeaf(self):
        """
        (GeneralTree) -> bool
        Returns true if this general tree is a leaf.
        """
        return self._degree == 0

    def getDegree(self):
        """
        (GeneralTree) -> int
        Returns the degree of this general tree node.
        """
        return self._degree

#{
    def getKey(self):
        """
        (GeneralTree) -> Object
        Returns the key in this general tree node.
        """
        return self._key

    def getSubtree(self, i):
        """
        (GeneralTree) -> Object
        Returns the specified subtree of this general tree node.
        """
        if i < 0 or i >= self._degree:
            raise IndexError
        ptr = self._list.head
        for j in xrange(i):
            ptr = ptr.next
        return ptr.datum

    def attachSubtree(self, t):
        """
        (GeneralTree, GeneralTree) -> None
        Attaches the given general tree as a subtree
        of this general tree node.
        """
        self._list.append(t)
        self._degree += 1

    def detachSubtree(self, t):
        """
        (GeneralTree, GeneralTree) -> GeneralTree
        Detaches and returns specified general tree
        from this general tree node.
        """
        self._list.extract(t)
        self._degree -= 1
        return t
#}>b

    def _compareTo(self, obj):
        """
        (GeneralTree, GeneralTree) -> int

        Compares this general tree with the given general tree.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "GeneralTree test program."
        print GeneralTree.main.__doc__
        gt = GeneralTree('A')
        gt.attachSubtree(GeneralTree('B'))
        gt.attachSubtree(GeneralTree('C'))
        gt.attachSubtree(GeneralTree('D'))
        gt.attachSubtree(GeneralTree('E'))
        Tree.test(gt)
        return 0

if __name__ == "__main__":
    sys.exit(GeneralTree.main(*sys.argv))
