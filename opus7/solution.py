#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: solution.py,v $
#   $Revision: 1.13 $
#
#   $Id: solution.py,v 1.13 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Solution class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.13 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.object import Object

#{
class Solution(Object):
    """
    Base class from which all solution space nodes are derived.
    """

    def __init__(self):
        """
        (Solution) -> None
        Constructor.
        """
        super(Solution, self).__init__()

    @abstractmethod
    def getIsFeasible(self): pass

    isFeasible = property(
        fget = lambda self: self.getIsFeasible())

    @abstractmethod
    def getIsComplete(self): pass

    isComplete = property(
        fget = lambda self: self.getIsComplete())

    @abstractmethod
    def getObjective(self): pass

    objective = property(
        fget = lambda self: self.getObjective())

    @abstractmethod
    def getBound(self): pass

    bound = property(
        fget = lambda self: self.getBound())

    @abstractmethod
    def getSuccessors(self): pass

    successors = property(
        fget = lambda self: self.getSuccessors())
#}>a
