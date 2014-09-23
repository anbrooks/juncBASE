#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: stackAsArray.py,v $
#   $Revision: 1.30 $
#
#   $Id: stackAsArray.py,v 1.30 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the StackAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.30 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.stack import Stack
from opus7.array import Array
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class StackAsArray(Stack):
    """
    Stack implemented using an array.
    """

#}@head

#{

    #...
#}@tail

#{
    def __init__(self, size = 0):
        """
        (StackAsArray [, int]) -> None
        Constructs a stack of the given size.
        """
        super(StackAsArray, self).__init__()
        self._array = Array(size)

    def purge(self):
        """
        (StackAsArray) -> None
        Purges this stack.
        """
        while self._count > 0:
            self._array[self._count] = None
            self._count -= 1
#}>a

#{
    def push(self, obj):
        """
        (StackAsArray, Object) -> None
        Pushes the given object on to this stack.
        """
        if self._count == len(self._array):
            raise ContainerFull
        self._array[self._count] = obj
        self._count += 1

    def pop(self):
        """
        (StackAsArray) -> Object
        Pops the top object off this stack.
        """
        if self._count == 0:
            raise ContainerEmpty
        self._count -= 1
        result = self._array[self._count]
        self._array[self._count] = None
        return result

    def getTop(self):
        """
        (StackAsArray) -> Object
        Returns the object at the top of this stack.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._array[self._count - 1]
#}>b

#{
    def accept(self, visitor):
        """
        (StackAsArray, Visitor) -> None
        Makes the given visitor visit all the objects in this stack.
        """
        assert isinstance(visitor, Visitor)
        for i in xrange(self._count):
            visitor.visit(self._array[i])
            if visitor.isDone:
                return
#}>c

#{
    class Iterator(Iterator):
        """
        Enumerates the elements of a StackAsArray.
        """

        def __init__(self, stack):
            """
            (StackAsArray.Iterator, Stack) -> None
            Constructs an iterator for the given stack.
            """
            super(StackAsArray.Iterator, self).__init__(stack)
            self._position = 0

        def next(self):
            """
            (StackAsArray.Iterator) -> Object
            Returns the next element.
            """
            if self._position >= self._container._count:
                raise StopIteration
	    obj = self._container._array[self._position]
            self._position = self._position + 1
            return obj

    def __iter__(self):
        """
        (StackAsArray) -> StackAsArray.Iterator
        Returns an iterator for this stack.
        """
        return self.Iterator(self)
#}>d
    
    def getIsFull(self):
        """
        (StackAsArray) -> bool
        Returns true of this stack is full.
        """
        return self._count == len(self._array)

    def _compareTo(self, obj):
        """
        (StackAsArray, StackAsArray) -> int

        Compares this stack with the given stack.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError
    
    @staticmethod
    def main(*argv):
        "StackAsArray test program."
        print StackAsArray.main.__doc__
        stack1 = StackAsArray(5)
        Stack.test(stack1)
        return 0

if __name__ == "__main__":
    sys.exit(StackAsArray.main(*sys.argv))
