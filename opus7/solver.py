#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: solver.py,v $
#   $Revision: 1.16 $
#
#   $Id: solver.py,v 1.16 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Solver class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.object import Object
from opus7.solution import Solution

#{
class Solver(Object):
    """
    Base class from which all problem solvers are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (Solver) -> None
        Constructor.
        """
        super(Solver, self).__init__()
        self._bestSolution = None
        self._bestObjective = sys.maxint

    @abstractmethod
    def search(self, initial): pass

    def solve(self, initial):
        """
        (Solver, Solution) -> Solution
        Solves a problem by searching the solution space
        starting from the given node.
        """
        assert isinstance(initial, Solution)
        self._bestSolution = None
        self._bestObjective = sys.maxint
        self.search(initial)
        return self._bestSolution

    def updateBest(self, solution):
        """
        (Solver, Solution) -> None
        Records the given solution if it is complete, feasible,
        and has a lower objective function value than the best
        solution seen so far.
        """
        if solution.isComplete and solution.isFeasible and \
                solution.objective < self._bestObjective:
            self._bestSolution = solution
            self._bestObjective = solution.objective
#}>a

    def _compareTo(self, obj):
        """
        (Solver, Solver) -> int

        Compares this solver with the given solver.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

