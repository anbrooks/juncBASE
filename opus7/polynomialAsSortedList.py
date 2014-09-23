#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: polynomialAsSortedList.py,v $
#   $Revision: 1.13 $
#
#   $Id: polynomialAsSortedList.py,v 1.13 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PolynomialAsSortedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.13 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.polynomial import Polynomial
from opus7.sortedListAsLinkedList import SortedListAsLinkedList
from opus7.visitor import Visitor
from copy import copy

#{
class PolynomialAsSortedList(Polynomial):
    """
    Polynomial implemented as a sorted list of terms.
    """

#}@head

#{

    # ...
#}@tail

    def addTerm(self, term):
        """
        (PolynomialAsSortedList, Polynomial.Term) -> None
        Adds the given term to this polynomial.
        """
        self._list.insert(term)

    def accept(self, visitor):
        """
        (PolynomialAsSortedList, Visitor) -> None
        Makes the given visitor visit all the terms in this polynomial.
        """
        assert isinstance(visitor, Visitor)
        self._list.accept(visitor)

    def find(self, term):
        """
        (PolynomialAsSortedList, Polynomial.Term) -> Polynomial.Term
        Finds a term in this polynomial that matches the given term.
        """
        return self._list.find(term)

    def withdraw(self, term):
        """
        (PolynomialAsSortedList, Polynomial.Term) -> Polynomial.Term
        Withdraws the given term from this polynomial.
        """
        self._list.withdraw(term)

#{
    def __init__(self):
        """
        (PolynomialAsSortedList) -> None
        Constructor.
        """
        super(PolynomialAsSortedList, self).__init__()
        self._list = SortedListAsLinkedList()

    def nextTerm(self, iter):
        """
        (PolynomialAsSortedList, Iterator) -> Polynomial.Term
        Returns the next term or None if there are no more terms.
        """
        try:
            return iter.next()
        except StopIteration:
            return None

    def __add__(self, poly):
        """
        (PolynomialAsSortedList, PolynomialAsSortedList)
            -> PolynomialAsSortedList.
        Returns the sum of this polynomial and the given polynomial.
        """
        result = PolynomialAsSortedList()
        p1 = iter(self._list)
        p2 = iter(poly._list)
        term1 = self.nextTerm(p1)
        term2 = self.nextTerm(p2)
        while term1 is not None and term2 is not None:
            if term1.exponent < term2.exponent:
                result.addTerm(copy(term1))
                term1 = self.nextTerm(p1)
            elif term1.exponent > term2.exponent:
                result.addTerm(copy(term2))
                term2 = self.nextTerm(p2)
            else:
                sum = term1 + term2
                if sum.coefficient != 0:
                    result.addTerm(sum)
                term1 = self.nextTerm(p1)
                term2 = self.nextTerm(p2)
        while term1 is not None:
            result.addTerm(copy(term1))
            term1 = self.nextTerm(p1)
        while term2 is not None:
            result.addTerm(copy(term2))
            term2 = self.nextTerm(p2)
        return result
#}>a
    
    def purge(self):
        """
        (PolynomialAsSortedList) -> None

        Purges this polynomial.
        """
        self._list.purge()

    def __iter__(self):
        raise NotImplementedError

    def _compareTo(self, obj):
        """
        (PolynomialAsSortedList, PolynomialAsSortedList) -> int

        Compares this polynomial with the given polynomial.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    def __str__(self):
        """
        (PolynomialAsSortedList) -> str
        Returns the string representation of this polynomial.
        """
        return str(self._list)
