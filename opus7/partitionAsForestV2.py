#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: partitionAsForestV2.py,v $
#   $Revision: 1.11 $
#
#   $Id: partitionAsForestV2.py,v 1.11 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PartitionaAsForestV2 class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.11 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.partition import Partition
from opus7.partitionAsForest import PartitionAsForest

#{
class PartitionAsForestV2(PartitionAsForest):
    """
    Partition (set of sets) implemented as a forest of trees.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, n):
        """
        (PartitionAsForestV2, int) -> None
        Constructs a partition with the given universe size.
        """
        super(PartitionAsForestV2, self).__init__(n)

#{
    def find(self, item):
        """
        (PartitionAsForestV2, item) -> Set
        Returns the element of this partition that contains the given item.
        (Collapsing find).
        """
        root = self._array[item]
        while root._parent is not None:
            root = root._parent
        ptr = self._array[item]
        while ptr._parent is not None:
            tmp = ptr._parent
            ptr._parent = root
            ptr = tmp
        return root
#}>a

#{
    def join(self, s, t):
        """
        (PartitionAsForestV2, Set, Set) -> None
        Joins the given elements of this partition.
        (Union by size).
        """
	assert s in self and s._parent is None and \
	    t in self and t._parent is None and s is not t
        if s._count > t._count:
            t._parent = s
            s._count += t._count
        else:
            s._parent = t
            t._count += s._count
        self._count -= 1
#}>b

    @staticmethod
    def main(*argv):
        "PartitionAsForestV2 test program."
        print PartitionAsForestV2.main.__doc__
        Partition.test(PartitionAsForestV2(5))
        return 0

if __name__ == "__main__":
    sys.exit(PartitionAsForestV2.main(*sys.argv))
