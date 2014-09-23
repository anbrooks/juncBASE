#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: partitionAsForest.py,v $
#   $Revision: 1.24 $
#
#   $Id: partitionAsForest.py,v 1.24 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the PartitionAsForest class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.24 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.partition import Partition
from opus7.tree import Tree
from opus7.set import Set
from opus7.array import Array
from opus7.visitor import Visitor

#{
class PartitionAsForest(Partition):
    """
    A partition (set of sets) implemented as a forest of trees.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, n):
        """
        (PartitionAsForest, int) -> None
        Constructs a partition with the given universe size.
        """
        super(PartitionAsForest, self).__init__(n)
        self._array = Array(self._universeSize)
        for item in xrange(self._universeSize):
            self._array[item] = self.PartitionTree(self, item)
        self._count = self._universeSize
#}>b

    def purge(self):
        """
        (PartitionAsForest) -> None
        Purges this partition.
        """
        for item in xrange(self._universeSize):
            self._array[item].purge()

#{
    def find(self, item):
        """
        (PartitionAsForest, int) -> Set
        Finds the set in this partition that contains the given element.
        """
        ptr = self._array[item]
        while ptr._parent is not None:
            ptr = ptr._parent
        return ptr
#}>c

#{
    def join(self, s, t):
        """
        (PartitionAsForest, Set, Set) -> None
        Joins the given sets in this partition.
        """
	assert s in self and s._parent is None and \
	    t in self and t._parent is None and s is not t
        t._parent = s
        self._count -= 1
#}>d

    def __contains__(self, obj):
        """
        (PartitionAsForest, Set) -> bool
        Returns true if the given set is an element of this partition.
        """
        return obj.isMemberOf(self)

    def accept(self, visitor):
        """
        (PartitionAsForest, Visitor) -> None

        Makes the given visitor visit the elements of this partition.
        """
        assert isinstance(visitor, Visitor)
        for i in xrange(self._universeSize):
            visitor.visit(self._array[i])
            if visitor.isDone:
                return

#{
    class PartitionTree(Set, Tree):
        """
        Represents a element of a partition.
        """

#}@head

#{

        # ...
#}@tail

#{
        def __init__(self, partition, item):
            """
            (PartitionAsForest.PartitionTree, PartitionAsForest, int) -> None
            Constructs an element of the given partition
            that contains the given item from the universal set.
            """
            super(PartitionAsForest.PartitionTree, self) \
                .__init__(partition._universeSize)
            self._partition = partition
            self._item = item
            self._parent = None
            self._rank = 0
            self._count = 1
#}>a

        def purge(self):
            """
            (PartitionAsForest.PartitionTree) -> None
            Purges element of the partition.
            """
            self._parent = None
            self._rank = 0
            self._count = 1

        def getCount(self):
            """
            (PartitionAsForest.PartitionTree) -> None
            Returns the number of items in this element of the partition.
            """
            return self._count
        
        def isMemberOf(self, partition):
            """
            (PartitionAsForest.PartitionTree, PartitionAsForest) -> None
            Returns true if this is an element of the given partition.
            """
            return self._partition is partition

        def getHeight(self):
            """
            (PartitionAsForest.PartitionTree) -> int
            Returns the height of this partition tree.
            """
            return self._rank

        def getKey(self):
            """
            (PartitionAsForest.PartitionTree) -> int
            Returns the item at the root of this partition tree.
            """
            return self._item

        def _compareTo(self, tree):
            """
            (PartitionAsForest.PartitionTree) -> int
            Compares the item at the root of this partition tree
            with the item at the root of the given partition tree.
            """
            assert isinstance(self, tree.__class__)
            return self._item - tree._item

        def __hash__(self):
            """
            (PartitionAsForest.PartitionTree) -> int
            Returns a hash of this partition tree.
            """
            return self._item

        def __str__(self):
            """
            (PartitionAsForest.PartitionTree) -> str
            Returns a string representation of this partition tree.
            """
            result = "PartitionTree {" + str(self._item)
            if self._parent is not None:
                result =  result + ", " + str(self._parent)
            return result + "}"

        def getIsEmpty(self):
            """
            (PartitionAsForest.PartitionTree) -> bool
            Returns false always.
            """
            return False

        def __iter__(self):
            raise NotImplementedError
        def __contains__(self):
            raise NotImplementedError
        def getIsLeaf(self):
            raise NotImplementedError
        def getDegree(self):
            raise NotImplementedError
        def getSubtree(self, i):
            raise NotImplementedError
        def insert(self, obj):
            raise NotImplementedError
        def withdraw(self, obj):
            raise NotImplementedError
        def __or__(self, tree):
            raise NotImplementedError
        def __and__(self, tree):
            raise NotImplementedError
        def __sub__(self, tree):
            raise NotImplementedError
        def __le__(self, tree):
            raise NotImplementedError
        def __eq__(self, tree):
            raise NotImplementedError

    def __iter__(self):
        raise NotImplementedError
    def insert(self, obj):
        raise NotImplementedError
    def withdraw(self, obj):
        raise NotImplementedError
    def __or__(self, partition):
        raise NotImplementedError
    def __and__(self, partition):
        raise NotImplementedError
    def __sub__(self, partition):
        raise NotImplementedError
    def __le__(self, partition):
        raise NotImplementedError
    def __eq__(self, partition):
        raise NotImplementedError

    def _compareTo(self, obj):
        """
        (PartitionAsForest, PartitionAsForest) -> int

        Compares this partition with the given partition.
        """
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "PartitionAsForest test program."
        print PartitionAsForest.main.__doc__
        Partition.test(PartitionAsForest(5))
        return 0

if __name__ == "__main__":
    sys.exit(PartitionAsForest.main(*sys.argv))
