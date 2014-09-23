#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: stackAsLinkedList.py,v $
#   $Revision: 1.32 $
#
#   $Id: stackAsLinkedList.py,v 1.32 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the StackAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.32 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.stack import Stack
from opus7.linkedList import LinkedList
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class StackAsLinkedList(Stack):
    """
    Stack implemented using a linked list.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (StackAsLinkedList)
        Constructs a stack.
        """
        super(StackAsLinkedList, self).__init__()
        self._list = LinkedList()

    def purge(self):
        """
        (StackAsLinkedList) -> None
        Purges this stack.
        """
        self._list.purge()
        self._count = 0
#}>a

#{
    def push(self, obj):
        """
        (StackAsLinkedList, Object) -> None
        Pushes the given object on to this stack.
        """
        self._list.prepend(obj)
        self._count += 1

    def pop(self):
        """
        (StackAsLinkedList) -> None
        Pops the top object off this stack.
        """
        if self._count == 0:
            raise ContainerEmpty
        result = self._list.first
        self._list.extract(result)
        self._count -= 1
        return result

    def getTop(self):
        """
        (StackAsLinkedList) -> None
        Returns the object at the top of this stack.
        """
        if self._count == 0:
            raise ContainerEmpty
        return self._list.first
#}>b

#{
    def accept(self, visitor):
        """
        (StackAsLinkedList, Visitor) -> None
        Makes the given visitor visit all the objects in this stack.
        """
        assert isinstance(visitor, Visitor)
        ptr = self._list.head
        while ptr is not None:
            visitor.visit(ptr.datum)
            if visitor.isDone:
                return
            ptr = ptr.next
#}>c

#{
    class Iterator(Iterator):
        """
        Enumerates the elements of a StackAsLinkedList.
        """

        def __init__(self, stack):
            """
            (StackAsLinkedList.Iterator, StackAsLinkedList) -> None
            Constructs an iterator for the given stack.
            """
            super(StackAsLinkedList.Iterator, self).__init__(
		stack)
            self._position = stack._list.head

        def next(self):
            """
            (StackAsLinkedList.Iterator) -> Object
            Returns the next element.
            """
            if self._position is None:
		raise StopIteration
	    element = self._position
	    self._position = self._position.next
            return element.datum

    def __iter__(self):
        """
        (StackAsLinkedList) -> StackAsLinkedList.Iterator
        Returns an iterator for this stack.
        """
        return self.Iterator(self)
#}>d
    
    def _compareTo(self, obj):
        """
        (StackAsLinkedList, StackAsLinkedList) -> int

        Compares this stack with the given stack.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "StackAsLinkedList test program."
        print StackAsLinkedList.main.__doc__
        stack2 = StackAsLinkedList()
        Stack.test(stack2)
        return 0

if __name__ == "__main__":
    sys.exit(StackAsLinkedList.main(*sys.argv))
