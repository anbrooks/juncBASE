#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: orderedListAsLinkedList.py,v $
#   $Revision: 1.31 $
#
#   $Id: orderedListAsLinkedList.py,v 1.31 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the OrderedListAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.31 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import exceptions
from opus7.orderedList import OrderedList
from opus7.linkedList import LinkedList
from opus7.cursor import Cursor
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.exception import *

#{
class OrderedListAsLinkedList(OrderedList):
    """
    Ordered list implemented using a linked list.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (OrderedListAsLinkedList) -> None
        Constructs an ordered list.
        """
        super(OrderedListAsLinkedList, self).__init__()
        self._linkedList = LinkedList()
#}>a

#{
    def insert(self, obj):
        """
        (OrderedListAsLinkedList, Object) -> None
        Inserts the given object at the end of this list.
        """
        self._linkedList.append(obj)
        self._count += 1

    def __getitem__(self, offset):
        """
        (OrderedListAsLinkedList, int) -> Object
        Returns the object in this list at the given offset.
        """
        if offset < 0 or offset >= self._count:
            raise IndexError
        ptr = self._linkedList.head
        i = 0
        while i < offset and ptr is not None:
            ptr = ptr.next
            i += 1
        return ptr.datum
#}>b

    def purge(self):
        """
        (OrderedListAsLinkedList) -> None
        Purges this ordered list.
        """
        self._linkedList = LinkedList()
        self._count = 0

    def accept(self, visitor):
        """
        (OrderedListAsLinkedList, Visitor) -> None
        Makes the given visitor visit all the objects in this ordered list.
        """
        assert isinstance(visitor, Visitor)
        ptr = self._linkedList.head
        while ptr is not None:
            visitor.visit(ptr.datum)
            if visitor.isDone:
                return
            ptr = ptr.next

#{
    def __contains__(self, obj):
        """
        (OrderedListAsLinkedList, Object) -> bool
        Returns true if the given object instance is in this ordered list.
        """
        ptr = self._linkedList.head
        while ptr is not None:
            if ptr.datum is obj:
                return True
            ptr = ptr.next
        return False

    def find(self, arg):
        """
        (OrderedListAsLinkedList, Object) -> Object
        Finds an object in this ordered list that equals the given object.
        """
        ptr = self._linkedList.head
        while ptr is not None:
            obj = ptr.datum
            if obj == arg:
                return obj
            ptr = ptr.next
        return None
#}>c

#{
    def withdraw(self, obj):
        """
        (OrderedListAsLinkedList, Object) -> None
        Withdraws the given object instance from this ordered list.
        """
        if self._count == 0:
            raise ContainerEmpty
        self._linkedList.extract(obj)
        self._count -= 1
#}>d

#{
    def findPosition(self, obj):
        """
        (OrderedListAsLinkedList, Object) -> OrderedListAsLinkedList.Cursor
        Finds the position of an object in this list
        that equals the given object and returns a cursor
        that refers to that object.
        """
        ptr = self._linkedList.head
        while ptr is not None:
            if ptr.datum == obj:
                break
            ptr = ptr.next
        return self.Cursor(self, ptr)
#}>e

#{
    class Cursor(Cursor):
        """
        A cursor that refers to an object in an ordered list.
        """

#}@head

#{

        # ...
#}@tail

#{
        def __init__(self, list, element):
            """
            (OrderedListAsLinkedList.Cursor, OrderedListAsLinkedList, 
                LinkedList.Element) -> None
            Constructs a cursor that refers to the object
           in the given element of the given list.
           """
	    super(OrderedListAsLinkedList.Cursor, self) \
		.__init__(list)
            self._element = element

        def getDatum(self):
            """
            (OrderedListAsLinkedList.Cursor) -> Object
            Returns the object to which this cursor refers.
            """
            return self._element.datum
#}>h

#{
        def insertAfter(self, obj):
            """
            (OrderedListAsLinkedList.Cursor, Object) -> None
            Inserts the given object into the list
            after the object to which this cursor refers.
            """
            self._element.insertAfter(obj)
            self._list._count += 1
#}>f

        def insertBefore(self, obj):
            """
            (OrderedListAsLinkedList.Cursor, Object) -> None
            Inserts the given object into the list
            before the object to which this cursor refers.
            """
            self._element.insertBefore(obj)
            self._list._count += 1

#{
        def withdraw(self):
            """
            (OrderedListAsLinkedList.Cursor) -> None
            Withdraws from the list the object to which this cursor refers.
            """
            self._list._linkedList.extract(self._element.datum)
            self._list._count -= 1
#}>g

    class Iterator(Iterator):
        """
        Enumerates the items in an ordered list.
        """

        def __init__(self, list):
            """
            (OrderedListAsLinkedList.Iterator, OrderedListAsLinkedList) -> None
            Constructs an iterator for the given list.
            """
            super(OrderedListAsLinkedList.Iterator, self).__init__(list)
            self._element = None

        def next(self):
            """
            (OrderedListAsLinkedList.Iterator) -> Object
            Returns the next element.
            """
            if self._element is None:
                self._element = self._container._linkedList.head
            else:
                self._element = self._element.next
            if self._element is None:
                raise StopIteration
            return self._element.datum

    def __iter__(self):
        """
        (OrderedListAsLinkedList) -> OrderedListAsLinkedList.Iterator
        Returns an iterator for this ordered list.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (OrderedListAsLinkedList, OrderedListAsLinkedList) -> int

        Compares this ordered list with the given ordered list.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "OrderedListAsLinkedList test program."
        print OrderedListAsLinkedList.main.__doc__
        list = OrderedListAsLinkedList()
        OrderedList.test(list)
        return 0

if __name__ == "__main__":
    sys.exit(OrderedListAsLinkedList.main(*sys.argv))
