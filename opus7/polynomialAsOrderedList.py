#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: polynomialAsOrderedList.py,v $
#   $Revision: 1.10 $
#
#   $Id: polynomialAsOrderedList.py,v 1.10 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PolynomialAsOrderedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.polynomial import Polynomial
from opus7.orderedListAsLinkedList import OrderedListAsLinkedList
from opus7.visitor import Visitor

#{
class PolynomialAsOrderedList(Polynomial):
    """
    Polynomial implemented as an ordered list of terms.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (PolynomialAsOrderedList) -> None
        Constructor.
        """
        self._list = OrderedListAsLinkedList()

    def addTerm(self, term):
        """
        (PolynomialAsOrderedList) -> None
        Adds the given term to this polynomial.
        """
        self._list.insert(term)

    def accept(self, visitor):
        """
        (PolynomialAsOrderedList, visitor) -> None
        Makes the given visitor visit all the terms in this polynomial.
        """
        assert isinstance(visitor, Visitor)
        self._list.accept(visitor)

    def find(self, term):
        """
        (PolynomialAsOrderedList, Polynomial.Term) -> Polynomial.Term
        Finds a term in this polynomial that matches the given term.
        """
        return self._list.find(term)

    def withdraw(self, term):
        """
        (PolynomialAsOrderedList, Polynomial.Term) -> None
        Withdraws the given term from this polynomial.
        """
        self._list.withdraw(term)
#}>a

    def purge(self):
        """
        (PolynomialAsOrderedList) -> None

        Purges this polynomial.
        """
        self._list.purge()

    def __iter__(self):
        raise NotImplementedError

    def __add__(self, polynomial):
        raise NotImplementedError

    def _compareTo(self, obj):
        """
        (PolynomialAsOrderedList, PolynomialAsOrderedList) -> int

        Compares this polynomial with the given polynomial.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    def __str__(self):
        """
        (PolynomialAsOrderedList) -> str
        Returns a string representation of this polynomial.
        """
        return str(self._list)
