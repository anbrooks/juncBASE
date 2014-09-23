#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: scalesBalancingProblem.py,v $
#   $Revision: 1.18 $
#
#   $Id: scalesBalancingProblem.py,v 1.18 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the ScalesBalancingProblem class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.18 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.solution import Solution
from opus7.solver import Solver
from opus7.array import Array
from opus7.iterator import Iterator
from copy import copy

class ScalesBalancingProblem(object):
    """
    Represents a scale-balancing problem.
    """

    def __init__(self, weight):
        """
        (ScalesBalancingProblem, Array) -> None
        Constructs a scales balancing problem with the given array of weights.
        """
        self._weight = weight
        self._numberOfWeights = len(weight)

    def solve(self, solver):
        """
        (ScalesBalancingProblem, Solver) -> Solution
        Solves this problem using the given solver.
        """
        assert isinstance(solver, Solver)
        return solver.solve(self.Node(self))

    class Node(Solution):
        """
        Represents a node in the solution space
        of this scales balancing problem.
        """

        def __init__(self, problem):
            """
            (ScalesBalancingProblem.Node, ScalesBalancingProblem) -> None
            Constructs the initial node in the solution space
            of this scales balancing problem.
            """
            self._problem = problem
            self._diff = 0
            self._unplacedTotal = 0
            self._numberPlaced = 0
            self._pan = Array(self._problem._numberOfWeights)

            for i in xrange(self._problem._numberOfWeights):
                self._unplacedTotal += self._problem._weight[i]

        def __copy__(self):
            """
            (ScalesBalancingProblem.Node) -> ScalesBalancingProblem.Node
            Returns a shallow copy of this node in the solution space
            of this scales balancing problem.
            """
            result = ScalesBalancingProblem.Node(self._problem)
            result._diff = self._diff
            result._numberPlaced = self._numberPlaced
            for i in xrange(self._problem._numberOfWeights):
                result._pan[i] = self._pan[i]
            result._unplacedTotal = self._unplacedTotal
            return result

        def getObjective(self):
            """
            (ScalesBalancingProblem.Node) -> int
            Returns the value of the objective function for this node
            in the solution space of this scales balancing problem.
            """
            return abs(self._diff)

        def getBound(self):
            """
            (ScalesBalancingProblem.Node) -> int
            Returns a bound on the objective function
            for all the nodes below this node (including this node)
            in the solution space of this scales balancing problem.
            """
            if abs(self._diff) > self._unplacedTotal:
                return abs(self._diff) - self._unplacedTotal
            else:
                return 0

        def getIsFeasible(self):
            """
            (ScalesBalancingProblem.Node) -> bool
            Returns true if this node is a feasible node
            in the solution space of this scales balancing problem.
            """
            return True

        def getIsComplete(self):
            """
            (ScalesBalancingProblem.Node) -> bool
            Returns true if this node is a complete node
            in the solution space of this scales balancing problem.
            A node is complete when all weights have been placed in a pan.
            """
            return self._numberPlaced == self._problem._numberOfWeights

        def placeNext(self, pan):
            """
            (ScalesBalancingProblem.Node, int) -> bool
            Modifies this node in the solution space of the scales
            balancing problem by placing the next unplaced weight
            into the specified pan.
            """
            self._pan[self._numberPlaced] = pan
            if pan == 0:
                self._diff += self._problem._weight[self._numberPlaced]
            else:
                self._diff -= self._problem._weight[self._numberPlaced]
            self._unplacedTotal -= self._problem._weight[self._numberPlaced]
            self._numberPlaced += 1

        def __str__(self):
            """
            (ScalesBalancingProblem.Node) -> str
            Returns a string representing this node
            in the solution space of this scales balancing problem.
            """
            s = ""
            comma = False
            for i in xrange(self._numberPlaced):
                if comma:
                    s = s + ", "
                s = s + str(self._pan[i])
                comma = True
            s = s + ", diff = " + str(self._diff)
            return s

        class SuccessorIterator(Iterator):
            """
            Enumerates the successor nodes of a given node
            in the solution space of this scales balancing problem.
            """

            def __init__(self, node):
                """
                (ScalesBalancingProblem.Node.SuccessorIterator,
                    ScalesBalancingProblem.Node) -> None
                Constructs an iterator that enumerates the successors
                of a given node in the solution space
                of this scales balancing problem.
                """
                super(ScalesBalancingProblem.Node.SuccessorIterator, self). \
                    __init__(node._problem)
                self._node = node
                self._pan = -1

            def next(self):
                """
                (ScalesBalancingProblem.Node.SuccessorIterator) ->
                    ScalesBalancingProblem.Node
                Returns the next successor of this node 
                in the solution space of this scales balancing problem.
                """
                self._pan += 1
                if self._pan == 2:
                    self._pan = -1
                    raise StopIteration
                result = copy(self._node)
                result.placeNext(self._pan)
                return result

        def getSuccessors(self):
            """
            (ScalesBalancingProblem.Node) ->
                ScalesBalancingProblem.Node.Iterator
            Returns an iterator that enumerates the successors of this node
            in the solution space of this scales balancing problem.
            """
            return self.SuccessorIterator(self)

        def _compareTo(self, obj):
            """
            (ScalesBalancingProblem.Node, ScalesBalancingProblem.Node) -> int

            Compares this noe with the given noe.
            """
            assert isinstance(self, obj.__class__)
            raise NotImplementedError
