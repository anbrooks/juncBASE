#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sortedListAsLinkedList.py,v $
#   $Revision: 1.20 $
#
#   $Id: sortedListAsLinkedList.py,v 1.20 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SortedListAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.20 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.sortedList import SortedList
from opus7.orderedListAsLinkedList import OrderedListAsLinkedList

#{
class SortedListAsLinkedList(
    OrderedListAsLinkedList, SortedList):
    """
    A sorted list implemented using a linked list.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (SortedListAsLinkedList) -> None
        Constructs a sorted list.
        """
        super(SortedListAsLinkedList, self).__init__()
#}>a

#{
    def insert(self, obj):
        """
        (SortedListAsLinkedList, Object) -> None
        Inserts the given object into this sorted list.
        """
        prevPtr = None
        ptr = self._linkedList.head
        while ptr is not None:
            if ptr.datum >= obj:
                break
            prevPtr = ptr
            ptr = ptr.next
        if prevPtr is None:
            self._linkedList.prepend(obj)
        else:
            prevPtr.insertAfter(obj)
        self._count += 1
#}>b

    def findElement(self, obj):
        """
        (SortedListAsLinkedList, Object) -> LinkedList.Element
        Finds the list element that contains an object
        that equals the given object.
        """
        ptr = self._linkedList.head
        while ptr is not None:
            if ptr.datum == obj:
                return ptr
            ptr = ptr.next
        return None

    class Cursor(OrderedListAsLinkedList.Cursor):
        """
        A cursor that refers to an object in a sorted list.
        """

        def __init__(self, list, element):
            """
            (SortedListAsLinkedList.Cursor, SortedListAsLinkedList,
                LinkedList.Element) -> None
            Constructs a cursor that refers to the object
            in the given element of the given list.
            """
            super(SortedListAsLinkedList.Cursor, self) \
                .__init_(list, element)

        def insertAfter(self, obj):
            """
            (SortedListAsLinkedList.Cursor, Object) -> None
            Not allowed in sorted lists.
            """
            raise TypeError

        def insertBefore(self, obj):
            """
            (SortedListAsLinkedList.Cursor, Object) -> None
            Not allowed in sorted lists.
            """
            raise TypeError

    def findPosition(self, obj):
        """
        (SortedListAsLinkedList, Object) -> SortedListAsLinkedList.Cursor
        Finds the position of an object in this sorted list
        that equals the given object and returns a cursor
        that refers to that object.
        """
        return self.Cursor(self, self.findElement(obj))
    
    @staticmethod
    def main(*argv):
        "SortedListAsLinkedList test program."
        print SortedListAsLinkedList.main.__doc__
        slist = SortedListAsLinkedList()
        SortedList.test(slist)
        return 0

if __name__ == "__main__":
    sys.exit(SortedListAsLinkedList.main(*sys.argv))
