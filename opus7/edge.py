#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: edge.py,v $
#   $Revision: 1.14 $
#
#   $Id: edge.py,v 1.14 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Edge class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.14 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.object import Object

#{
class Edge(Object):
    """
    Base class from which all graph edge classes are derived.
    """

    def __init__(self):
        """
        (Edge) -> None
        Constructor.
        """
        super(Edge, self).__init__()

    @abstractmethod
    def getV0(self): pass

    v0 = property(
        fget = lambda self: self.getV0())

    @abstractmethod
    def getV1(self): pass

    v1 = property(
        fget = lambda self: self.getV1())

    @abstractmethod
    def getWeight(self): pass

    weight = property(
        fget = lambda self: self.getWeight())

    @abstractmethod
    def getIsDirected(self): pass

    isDirected = property(
        fget = lambda self: self.getIsDirected())

    @abstractmethod
    def mateOf(self, vertex): pass
#}>a
