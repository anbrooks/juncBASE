#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: zeroOneKnapsackProblem.py,v $
#   $Revision: 1.17 $
#
#   $Id: zeroOneKnapsackProblem.py,v 1.17 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the ZeroOneKnapsackProblem class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.17 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from copy import copy
from opus7.solution import Solution
from opus7.solver import Solver
from opus7.array import Array
from opus7.iterator import Iterator

class ZeroOneKnapsackProblem(object):
    """
    Represents a zero-one knapsack problem.
    """

    def __init__(self, weight, profit, capacity):
        """
        (ZeroOneKnapsackProblem, Array, Array, int) -> None
        Constructs a zero-one knapsack problem with the given
        array of weights, profits, and capacity.
        """
        self._numberOfItems = len(weight)
        self._weight = weight
        self._profit = profit
        self._capacity = capacity

    def solve(self, solver):
        """
        (ZeroOneKnapsackProblem, Solver) -> Solution
        Solves this problem using the given solver.
        """
        assert isinstance(solver, Solver)
        return solver.solve(self.Node(self))

    class Node(Solution):
        """
        Represents a node in the solution space
        of this zero-one knapsack problem.
        """

        def __init__(self, problem):
            """
            (ZeroOneKnapsackProblem.Node, ZeroOneKnapsackProblem) -> None
            Constructs the initial node in the solution space
            of this zero-one knapsack problem.
            """
            self._problem = problem
            self._totalWeight = 0
            self._totalProfit = 0
            self._unplacedProfit = 0
            self._numberPlaced = 0
            self._x = Array(problem._numberOfItems)
            for i in xrange(self._problem._numberOfItems):
                self._unplacedProfit += self._problem._profit[i]

        def __copy__(self):
            """
            (ZeroOneKnapsackProblem.Node) -> ZeroOneKnapsackProblem
            Returns a shallow copy of this node in the solution space
            of this zero-one knapsack problem.
            """
            result = ZeroOneKnapsackProblem.Node(self._problem)
            result._totalWeight = self._totalWeight
            result._totalProfit = self._totalProfit
            result._numberPlaced = self._numberPlaced
            for i in xrange(self._problem._numberOfItems):
                result._x[i] = self._x[i]
            result._unplacedProfit = self._unplacedProfit
            return result

        def getObjective(self):
            """
            (ZeroOneKnapsackProblem.Node) -> int
            Returns the value of the objective function for this node
            in the solution space of this zero-one knapsack problem.
            """
            return -self._totalProfit

        def getBound(self):
            """
            (ZeroOneKnapsackProblem.Node) -> int
            Returns a bound on the value of the objective function
            for all the nodes below this node (including this node)
            in the solution space of this zero-one knapsack problem.
            """
            return -(self._totalProfit + self._unplacedProfit)

        def getIsFeasible(self):
            """
            (ZeroOneKnapsackProblem.Node) -> bool
            Returns true if this node is a feasible node
            in the solution space of this zero-one knapsack problem.
            """
            return self._totalWeight <= self._problem._capacity

        def getIsComplete(self):
            """
            (ZeroOneKnapsackProblem.Node) -> bool
            Returns true if this node is a feasible node
            in the solution space of this zero-one knapsack problem.
            A node is complete when no more items can be placed in the knapsack.
            """
            return self._numberPlaced == self._problem._numberOfItems

        def placeNext(self, value):
            """
            (ZeroOneKnapsackProblem.Node, bool) -> None
            Modifies this node in the solution space of
            the zero-one knapsack problem by adding the next unplaced
            item into the knapsack if value is true.
            """
            self._x[self._numberPlaced] = value
            if value == 1:
                self._totalWeight += self._problem._weight[self._numberPlaced]
                self._totalProfit += self._problem._profit[self._numberPlaced]
                self._unplacedProfit -= \
                    self._problem._profit[self._numberPlaced]
            self._numberPlaced += 1

        def __str__(self):
            """
            (ZeroOneKnapsackProblem.Node, bool) -> str
            Returns a string representation of this node
            in the solution space of this zero-one knapsack problem.
            """
            comma = False
            s = ""
            for i in xrange(self._numberPlaced):
                if comma:
                    s = s + ", "
                s = s + str(self._x[i])
                comma = True
            s = s + ", total weight = " + str(self._totalWeight)
            s = s + ", total profit = " + str(self._totalProfit)
            return s

        class SuccessorIterator(Iterator):
            """
            Enumerates the successor nodes of this node
            in the solution space of this zero-one knapsack problem.
            """

            def __init__(self, node):
                """
                (ZeroOneKnapsackProblem.Node.SuccessorIterator,
                    ZeroOneKnapsackProblem.Node) -> None
                Constructs an iterator that enumerates the successors
                of a given node in the solution space
                of this zero-one knapsack problem.
                """
                super(ZeroOneKnapsackProblem.Node.SuccessorIterator, self) \
                    .__init__(node._problem)
                self._node = node
                self._x = -1

            def next(self):
                """
                (ZeroOneKnapsackProblem.Node.SuccessorIterator) ->
                    ZeroOneKnapsackProblem.Node
                Returns the next succeessor of this node
                in the solution space of this zero-one knapsack problem.
                """
                self._x += 1
                if self._x == 2:
                    self._x = -1
                    raise StopIteration
                result = copy(self._node)
                result.placeNext(self._x)
                return result

        def getSuccessors(self):
            """
            (ZeroOneKnapsackProblem.Node)
                -> ZeroOneKnapsackProblem.Node.SuccessorIterator
            Returns an iterator that enumerates the successors of this node
            in the solution space of this zero-one knapsack problem.
            """
            return self.SuccessorIterator(self)

        def _compareTo(self, obj):
            """
            (ZeroOneKnapsackProblem.Node, ZeroOneKnapsackProblem.Node) -> int

            Compares this node with the given node.
            """
            assert isinstance(self, obj.__class__)
            raise NotImplementedError
