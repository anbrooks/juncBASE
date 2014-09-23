#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: multisetAsLinkedList.py,v $
#   $Revision: 1.23 $
#
#   $Id: multisetAsLinkedList.py,v 1.23 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the MultisetAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.23 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.multiset import Multiset
from opus7.linkedList import LinkedList
from opus7.iterator import Iterator
from opus7.visitor import Visitor

#{
class MultisetAsLinkedList(Multiset):
    """
    Multiset implemented using a linked list of elements.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, n):
        """
        (MultisetAsLinkedList, int) -> None
        Constructs a multiset with the given universe size.
        """
        super(MultisetAsLinkedList, self).__init__(n)
        self._list = LinkedList()
#}>a

    def insert(self, item):
        """
        (MultisetAsLinkedList, Object) -> None
        Inserts the given element into this multiset.
        """
        ptr = self._list.head
        prevPtr = None
        while ptr is not None:
            if ptr.datum >= item:
                break
            prevPtr = ptr
            ptr = ptr.next
        if prevPtr is None:
            self._list.prepend(item)
        else:
            prevPtr.insertAfter(item)

    def withdraw(self, item):
        """
        (MultisetAsLinkedList, Object) -> None
        Withdraws the given element from this multiset.
        """
        ptr = self._list.head
        while ptr is not None:
            if ptr.datum == item:
                list.extract(ptr)
                return
            ptr = ptr.next

    def __contains__(self, item):
        """
        (MultisetAsLinkedList, Object) -> bool
        Returns true if the given elements is in this multiset.
        """
        ptr = self._list.head
        while ptr is not None:
            if ptr.datum == item:
                return True
            ptr = ptr.next
        return False

    def purge(self):
        """
        (MultisetAsLinkedList) -> None
        Purges this multiset.
        """
        self._list = LinkedList()

    def getCount(self):
        """
        (MultisetAsLinkedList) -> int
        Returns the number of elements in this multiset.
        """
        result = 0
        ptr = self._list.head
        while ptr is not None:
            result += 1
            ptr = ptr.next
        return result

    def accept(self, visitor):
        """
        (MultisetAsLinkedList, Visitor) -> None
        Makes the given visitor visit all the elements in this multiset.
        """
        assert isinstance(visitor, Visitor)
        ptr = self._list.head
        while ptr is not None:
            visitor.visit(ptr.datum)
            if visitor.isDone:
                return
            ptr = ptr.next

#{
    def __or__(self, set):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> MultisetAsLinkedList
        Returns the union of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsLinkedList)
	assert self._universeSize == set._universeSize
        result = MultisetAsLinkedList(self._universeSize)
        p = self._list.head
        q = set._list.head
        while p is not None and q is not None:
            if p.datum <= q.datum:
                result._list.append(p.datum)
                p = p.next
            else:
                result._list.append(q.datum)
                q = q.next
        while p is not None:
            result._list.append(p.datum)
            p = p.next
        while q is not None:
            result._list.append(q.datum)
            q = q.next
        return result
#}>b

#{
    def __and__(self, set):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> MultisetAsLinkedList
        Returns the intersection of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsLinkedList)
	assert self._universeSize == set._universeSize
        result = MultisetAsLinkedList(self._universeSize)
        p = self._list.head
        q = set._list.head
        while p is not None and q is not None:
            diff = p.datum - q.datum
            if diff == 0:
                result._list.append(p.datum)
            if diff <= 0:
                p = p.next
            if diff >= 0:
                q = q.next
        return result
#}>c

    def __sub__(self, set):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> MultisetAsLinkedList
        Returns the difference of this multiset and the given multiset.
        """
	assert isinstance(set, MultisetAsLinkedList)
	assert self._universeSize == set._universeSize
        result = MultisetAsLinkedList(self._universeSize)
        p = self._list.head
        q = set._list.head
        while p is not None and q is not None:
            diff = p.datum - q.datum
            if diff < 0:
                result._list.append(p.datum)
            if diff <= 0:
                p = p.next
            if diff >= 0:
                q = q.next
        while p is not None:
            result._list.append(p.datum)
            p = p.next
        return result

    def __le__(self, set):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> bool
        Returns true if this multiset is a proper subset of the given multiset.
        """
	assert isinstance(set, MultisetAsLinkedList)
	assert self._universeSize == set._universeSize
        p = self_list.head
        q = set._list.head
        while p is not None and q is not None:
            diff = p.datum - q.datum
            if diff == 0:
                p = p.next
                q = q.next
            elif diff > 0:
                q = q.next
            else:
                return False
        if p is not None:
            return False
        else:
            return True

    def __eq__(self, set):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> bool
        Returns true if this multiset is equal to the given multiset.
        """
	assert isinstance(set, MultisetAsLinkedList)
	assert self._universeSize == set._universeSize
        p = self._list.head
        q = set._list.head
        while p is not None and q is not None:
            if p.datum != q.datum:
                return False
            p = p.next
            q = q.next
        if p is not None or q is not None:
            return False
        else:
            return True

    class Iterator(Iterator):
        """
        Enumerates the elements of a MultisetAsLinkedList.
        """

        def __init__(self, multiset):
            """
            (MultisetAsLinkedList.Iterator, MultisetAsLinkedList) -> None
            Constructs an interator for the given multiset.
            """
            super(MultisetAsLinkedList.Iterator, self).__init__(multiset)
            self._ptr = None

        def next(self):
            """
            (MultisetAsLinkedList.Iterator) -> Object
            Returns the next element in the multiset.
            """
            if self._ptr is None:
                self._ptr = self._container._list.head
            else:
                self._ptr = self._ptr.next
            if self._ptr == None:
                raise StopIteration
            return self._ptr.datum

    def __iter__(self):
        """
        (MultisetAsLinkedList) -> MultisetAsLinkedList.Iterator
        Returns an iterator for this multiset.
        """
        return self.Iterator(self)

    def _compareTo(self, obj):
        """
        (MultisetAsLinkedList, MultisetAsLinkedList) -> int

        Compares this multiset with the given multiset.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError


    @staticmethod
    def main(*argv):
        "MultisetAsLinkedList test program."
        print MultisetAsLinkedList.main.__doc__
        Multiset.test(MultisetAsLinkedList(32), \
            MultisetAsLinkedList(32), MultisetAsLinkedList(32))
        return 0

if __name__ == "__main__":
    sys.exit(MultisetAsLinkedList.main(*sys.argv))
