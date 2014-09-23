#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: container.py,v $
#   $Revision: 1.41 $
#
#   $Id: container.py,v 1.41 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Container class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.41 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.object import Object
from opus7.visitor import Visitor

#{
class Container(Object):
    """
    Base class from which all containers are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (Container) -> None
        Constructs this container.
        """
        super(Container, self).__init__()
        self._count = 0

    @abstractmethod
    def purge(self): pass

    @abstractmethod
    def __iter__(self): pass
#}>a

#{
    def getCount(self):
        """
        (Container) -> int
        Returns the number of items in this container.
        """
        return self._count

    count = property(
        fget = lambda self: self.getCount())

    def getIsEmpty(self):
        """
        (Container) -> bool
        Returns true if this container is empty.
        """
        return self.count == 0

    isEmpty = property(
        fget = lambda self: self.getIsEmpty())

    def getIsFull(self):
        """
        (Container) -> bool
        Returns true if this container is full.
        """
        return False

    isFull = property(
        fget = lambda self: self.getIsFull())
#}>b

#{
    def accept(self, visitor):
        """
        (Container, Visitor) -> None
        Makes the given visitor visit all the items in this container.
        """
        assert isinstance(visitor, Visitor)
        for obj in self:
            visitor.visit(obj)
#}>c

    def elements(self):
        """
        (Container) -> Object

        Generator that yields the objects in this container.
        """
        for obj in self:
            yield obj

#{
    class StrVisitor(Visitor):
        """
        Visitor that accumulates visited items in a string.
        """

        def __init__(self):
            """
            (Container.StrVisitor) -> None
            Constructor.
            """
            self._string = ""
            self._comma = False

        def visit(self, object):
            """
            (Container.StrVisitor, Object) -> None
            Appends the given object to the string of visited items.
            """
            if self._comma:
                self._string = self._string + ", "
            self._string = self._string + str(object)
            self._comma = True

        def __str__(self):
            """
            (Container.StrVisitor) -> string
            Returns the string of visited items.
            """
            return self._string

    def __str__(self):
        """
        (Container) -> string
        Returns a string representation of this container.
        """
        visitor = Container.StrVisitor()
        self.accept(visitor)
        return "%s {%s}" % (self.__class__.__name__,str(visitor))
#}>d

#{
    def __hash__(self):
        """
        (Container) -> int
        Returns the hash of this container.
        """
	result = hash(self.__class__)
	for obj in self:
	    result = (result + hash(obj)) & sys.maxint
	return result
#}>e

    @staticmethod
    def main(*argv):
        "Container test program."
        print Container.main.__doc__

        class Dummy(Container):
            def _compareTo(self, obj): pass
            def purge(self): pass
            def __iter__(self): pass
        container = Dummy()
        print container._count
        return 0

if __name__ == "__main__":
    sys.exit(Container.main(*sys.argv))
