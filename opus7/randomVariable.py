#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: randomVariable.py,v $
#   $Revision: 1.17 $
#
#   $Id: randomVariable.py,v 1.17 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the RandomVariable class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.17 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.object import Object

#{
class RandomVariable(Object):
    """
    Base class from which all random variables are derived.
    """

    def __init__(self):
        """
        (RandomVariable) -> None
        Constructor.
        """
        super(RandomVariable, self).__init__()

    @abstractmethod
    def getNext(self): pass

    next = property(
        fget = lambda self: self.getNext())
#}>a

    def _compareTo(self, obj):
        """
        (RandomVariable, RandomVariable) -> int

        Compares this random variable with the given random variable.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

