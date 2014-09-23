#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: vertex.py,v $
#   $Revision: 1.14 $
#
#   $Id: vertex.py,v 1.14 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the Vertex class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.14 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.object import Object

#{
class Vertex(Object):
    """
    Base class from which all graph vertex classes are derived.
    """

    def __init__(self):
        """
        (Vertex) -> None
        Constructor.
        """
        super(Vertex, self).__init__()

    @abstractmethod
    def getNumber(self): pass

    number = property(
        fget = lambda self: self.getNumber())

    @abstractmethod
    def getWeight(self): pass

    weight = property(
        fget = lambda self: self.getWeight())

    @abstractmethod
    def getIncidentEdges(self): pass

    incidentEdges = property(
        fget = lambda self: self.getIncidentEdges())

    @abstractmethod
    def getEmanatingEdges(self): pass

    emanatingEdges = property(
        fget = lambda self: self.getEmanatingEdges())

    @abstractmethod
    def getPredecessors(self): pass

    predecessors = property(
        fget = lambda self: self.getPredecessors())

    @abstractmethod
    def getSuccessors(self): pass

    successors = property(
        fget = lambda self: self.getSuccessors())
#}>a
