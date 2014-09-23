#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: multiset.py,v $
#   $Revision: 1.8 $
#
#   $Id: multiset.py,v 1.8 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Multiset class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.8 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.object import Object
from opus7.set import Set

#{
class Multiset(Set):
    """
    Base class from which all multiset classes are derived.
    """

    def __init__(self, universeSize):
        """
        (Multiset, int) -> None
        Constructs a multiset with the given universe size.
        """
        super(Multiset, self).__init__(universeSize)
#}>a

    @staticmethod
    def test(s1, s2, s3):
        "Multiset test program."
        for i in xrange(0, 4):
            s1.insert(i)
        for i in xrange(2, 6):
            s2.insert(i)
        for i in xrange(0, 4, 2):
            s3.insert(i)
        print s1
        print s2
        print s3
        print s1 | s2 # union
        print s1 & s3 # intersection
        print s1 - s3 # difference
        for i in s3:
            print i,
        print
