#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: polynomial.py,v $
#   $Revision: 1.22 $
#
#   $Id: polynomial.py,v 1.22 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Polynomial class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.22 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.object import Object
from opus7.container import Container
from opus7.visitor import Visitor

#{
class Polynomial(Container):
    """
    Base class from which all polynomial classes are derived."
    """

#}@head

#{

    # ...
#}@tail


#{
    def __init__(self):
        """
        (Polynomail) -> None
        Constructor.
        """
        super(Polynomial, self).__init__()

    @abstractmethod
    def addTerm(self, term): pass

    @abstractmethod
    def differentiate(self): pass

    @abstractmethod
    def __add__(self, polynomial): pass
#}>a

#{
    class DifferentiatingVisitor(Visitor):
        """
        Visitor that differentiates the terms it visits.
        """

        def __init__(self):
            """
            (Polynomial.DifferentiatingVisitor) -> None
            Constructor.
            """
            super(Polynomial.DifferentiatingVisitor, self) \
                .__init__()

        def visit(self, term):
            """
            (Polynomial.DifferentiatingVisitor, Polynomial.Term)
                -> None
            Differentiates the given term.
            """
	    term.differentiate()

    def differentiate(self):
        """
        (Polynomial) -> None
        Differentiates this polynomial.
        """
        self.accept(self.DifferentiatingVisitor())
        zeroTerm = self.find(self.Term(0, 0))
        if zeroTerm is not None:
            self.withdraw(zeroTerm)
#}>d

#{
    class Term(Object):
        """
        Represents a term in a polynomial.
        """

#}@head

#{

        # ...
#}@tail

#{
        def __init__(self, coefficient, exponent):
            """
            (Polynomial.Term, double, int) -> None
            Constructs a term with the given coefficient and exponent.
            """
            self._coefficient = coefficient
            self._exponent = exponent

        def _compareTo(self, term):
            """
            (Polynomial.Term, Polynomial.Term) -> int
            Compares this term with the given term.
            """
            assert isinstance(self, term.__class__)
            if self._exponent == term._exponent:
		return cmp(self._coefficient, term._coefficient)
            else:
		return cmp(self._exponent, term._exponent)

        def differentiate(self):
            """
            (Polynomial.Term) -> None
            Differentiates this term.
            """
            if self._exponent > 0:
                self._coefficient *= self._exponent
                self._exponent -= 1
            else:
                self._coefficient = 0
#}>b

#{
        def __copy__(self):
            """
            (Polynomial.Term) -> Polynomial.Term
            Returns a shallow copy of this term.
            """
            return Polynomial.Term(
                self._coefficient, self._exponent)

        def getCoefficient(self):
            return self._coefficient

        coefficient = property(
            fget = lambda self: self.getCoefficient())

        def getExponent(self):
            return self._exponent

        exponent = property(
            fget = lambda self: self.getExponent())

        def __add__(self, term):
            """
            (Polynomial.Term, Polynomial.Term) -> Polynomial.Term
            Returns the sum of this term and the given term.
            """
	    assert self._exponent == term._exponent
            return Polynomial.Term(
                self._coefficient + term._coefficient,
                self._exponent)
#}>c

        def __str__(self):
            """
            (Polynomial.Term) -> str
            Returns a string representation of this term.
            """
            return "%gx^%g" % (self._coefficient, self._exponent)
