#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: partition.py,v $
#   $Revision: 1.13 $
#
#   $Id: partition.py,v 1.13 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Partition class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.13 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.set import Set

#{
class Partition(Set):
    """
    Base class from which all partitions are derived.
    """

    def __init__(self, n):
        super(Partition, self).__init__(n)

    @abstractmethod
    def find(self, obj): pass

    @abstractmethod
    def join(self, s, t): pass
#}>a

    @staticmethod
    def test(p):
        "Partition test program."
        print Partition.test.__doc__
        print p
        s2 = p.find(2)
        s4 = p.find(4)
        p.join(s2, s4)
        s3 = p.find(3)
        s4b = p.find(4)
        p.join(s3, s4b)
        print p
