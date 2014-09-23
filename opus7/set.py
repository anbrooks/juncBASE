#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: set.py,v $
#   $Revision: 1.20 $
#
#   $Id: set.py,v 1.20 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Set class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.20 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.searchableContainer import SearchableContainer

#{
class Set(SearchableContainer):
    """
    Base class from which all Set classes are derived.
    """

    def __init__(self, universeSize):
        """
        (Set, int) -> None
        Constructs a set with the given universe size.
        """
        self._universeSize = universeSize

    def getUniverseSize(self):
        return self._universeSize

    universeSize = property(
        fget = lambda self: self.getUniverseSize())

    @abstractmethod
    def __or__(self, set): pass

    @abstractmethod
    def __and__(self, set): pass

    @abstractmethod
    def __sub__(self, set): pass

    @abstractmethod
    def __eq__(self, set): pass

    @abstractmethod
    def __le__(self, set): pass
#}>a

    def find(self, i):
        """
        (Set, int) -> int

        Returns the given integer if it is in this set.
        Returns None otherwise.
        """
        if self.isMember(i):
            return i
        else:
            return None

    @staticmethod
    def test(s1, s2, s3):
        "Set test program."
        print Set.test.__doc__
        for i in xrange(0, 4, 1):
            s1.insert(i)
        for i in xrange(2, 6, 1):
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
