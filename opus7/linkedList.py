#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: linkedList.py,v $
#   $Revision: 1.25 $
#
#   $Id: linkedList.py,v 1.25 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the LinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.25 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.exception import *

#{
class LinkedList(object):
    """
    Linked list class.
    """

#}@head

#{

    # ...
#}@tail

#{
    class Element(object):
        """
        An element of a linked list.
        """

        def __init__(self, list, datum, next):
            """
            (LinkedList.Element, LinkedList, Object, LinkedList.Element) -> None
            Constructs a list element with the given values.
            """
            self._list = list
            self._datum = datum
            self._next = next

        def getDatum(self):
            """
            (LinkedList.Element) -> Object
            Returns the datum in this list element.
            """
            return self._datum

        datum = property(
            fget = lambda self: self.getDatum())

        def getNext(self):
            """
            (LinkedList.Element) -> LinkedList.Element
            Returns the next list element
            """
            return self._next

        next = property(
            fget = lambda self: self.getNext())
#}>a

#{
#!    class Element(object):

        def insertAfter(self, item):
            """
            (LinkedList.Element, Object) -> None
            Inserts the given item after this list element.
            """
            self._next = LinkedList.Element(
                self._list, item, self._next)
            if self._list._tail is self:
                self._list._tail = self._next

        def insertBefore(self, item):
            """
            (LinkedList.Element, Object) -> None
            Inserts the given item before this list element.
            """
            tmp = LinkedList.Element(self._list, item, self)
            if self is self._list._head:
                self._list._head = tmp
            else:
                prevPtr = self._list._head
                while prevPtr is not None \
                        and prevPtr._next is not self:
                    prevPtr = prevPtr._next
                prevPtr._next = tmp
#}>j

        def extract(self):
            """
            (LinkedList.Element) -> None
            Extracts this list element from the list.
            """
            prevPtr = None
            if self._list._head is self:
                self._list._head = next
            else:
                prevPtr = list._head
                while prevPtr is not None \
                        and prevPtr._next is not self:
                    prevPtr = prevPtr._next
                if prevPtr is None:
                    raise InternalError
                prevPtr._next = next
            if self._list._tail is self:
                self._list._tail = prevPtr

#{
    def __init__(self):
        """
        (LinkedList) -> None
        Constructs an empty linked list.
        """
        self._head = None
        self._tail = None
#}>b

#{
    def purge(self):
        """
        (LinkedList) -> None
        Purges this linked list.
        """
        self._head = None
        self._tail = None
#}>c

#{
    def getHead(self):
        """
        (LinkedList) -> LinkedList.Element
        Returns the list element at the head of this list.
        """
        return self._head

    head = property(
        fget = lambda self: self.getHead())
    
    def getTail(self):
        """
        (LinkedList) -> LinkedList.Element
        Returns the list element at the tail of this list.
        """
        return self._tail

    tail = property(
        fget = lambda self: self.getTail())
    
    def getIsEmpty(self):
        """
        (LinkedList) -> bool
        Returns true if this list is empty.
        """
        return self._head is None

    isEmpty = property(
        fget = lambda self: self.getIsEmpty())
#}>d

#{
    def getFirst(self):
        """
        (LinkedList) -> Object
        Returns the first item in this list.
        """
        if self._head is None:
            raise ContainerEmpty
        return self._head._datum

    first = property(
        fget = lambda self: self.getFirst())

    def getLast(self):
        """
        (LinkedList) -> Object
        Returns the last item in this list.
        """
        if self._tail is None:
            raise ContainerEmpty
        return self._tail._datum

    last = property(
        fget = lambda self: self.getLast())
#}>e

#{
    def prepend(self, item):
        """
        (LinkedList, Object) -> None
        Prepends the given item to this list.
        """
        tmp = self.Element(self, item, self._head)
        if self._head is None:
            self._tail = tmp
        self._head = tmp
#}>f

#{
    def append(self, item):
        """
        (LinkedList, Object) -> None
        Appends the given item to this list.
        """
        tmp = self.Element(self, item, None)
        if self._head is None:
            self._head = tmp
        else:
            self._tail._next = tmp
        self._tail = tmp
#}>g

#{
    def __copy__(self):
        """
        (LinkedList) -> LinkedList
        Returns a shallow copy of this linked list.
        """
        result = LinkedList()
        ptr = list._head
        while ptr is not None:
            result.append(ptr._datum)
            ptr = ptr._next
        return result
#}>h

#{
    def extract(self, item):
        """
        (LinkedList, Object) -> None
        Extracts the given item from this list.
        """
        ptr = self._head
        prevPtr = None
        while ptr is not None and ptr._datum is not item:
            prevPtr = ptr
            ptr = ptr._next
        if ptr is None:
            raise KeyError
        if ptr == self._head:
            self._head = ptr._next
        else:
            prevPtr._next = ptr._next
        if ptr == self._tail:
            self._tail = prevPtr
#}>i

    def __str__(self):
        """
        (LinkedList) -> string
        Returns a string representation of this list.
        """
        string = "LinkedList {"
        ptr = self._head
        while ptr is not None:
            string = string + str(ptr._datum)
            if ptr._next is not None:
                string = string + ", "
            ptr = ptr._next
        string = string + "}"
        return string

    @staticmethod
    def main(*argv):
        "LinkedList test program."
        print LinkedList.main.__doc__
        l1 = LinkedList()
        l1.append(57)
        l1.append("hello")
        l1.append(None)
        print l1
        print "isEmpty returns %s" % (l1.isEmpty)
        return 0

if __name__ == "__main__":
    sys.exit(LinkedList.main(*sys.argv))
