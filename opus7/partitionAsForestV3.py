#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: partitionAsForestV3.py,v $
#   $Revision: 1.10 $
#
#   $Id: partitionAsForestV3.py,v 1.10 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PartitionAsForestV2 class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.partition import Partition
from opus7.partitionAsForestV2 import PartitionAsForestV2

#{
class PartitionAsForestV3(PartitionAsForestV2):
    """
    Partition (set of sets) implemented as a forest of trees.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, n):
        """
        (PartitionAsForestV3, int) -> None
        Constructs a partition with the given universe size.
        """
        super(PartitionAsForestV3, self).__init__(n)

#{
    def join(self, s, t):
        """
        (PartitionAsForestV3, Set, Set) -> None
        Joins the given sets of this partition.
        (Union by rank).
        """
	assert s in self and s._parent is None and \
	    t in self and t._parent is None and s is not t
        if s._rank > t._rank:
            t._parent = s
        else:
            s._parent = t
            if s._rank == t._rank:
                t._rank += 1
        self._count -= 1
#}>a

    @staticmethod
    def main(*argv):
        "PartitionAsForestV3 test program."
        print PartitionAsForestV3.main.__doc__
        Partition.test(PartitionAsForestV3(5))
        return 0

if __name__ == "__main__":
    sys.exit(PartitionAsForestV3.main(*sys.argv))
