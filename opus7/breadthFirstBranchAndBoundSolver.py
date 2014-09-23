#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: breadthFirstBranchAndBoundSolver.py,v $
#   $Revision: 1.4 $
#
#   $Id: breadthFirstBranchAndBoundSolver.py,v 1.4 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the BreadthFirstBrandAndBoundSolver class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.4 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.solver import Solver
from opus7.queueAsLinkedList import QueueAsLinkedList

class BreadthFirstBranchAndBoundSolver(Solver):
    """
    Breadth-first brand-and-bound solver.
    """

    def __init__(self):
        """
        (BreadthFirstBranchAndBoundSolver) -> None
        Constructor.
        """
        super(BreadthFirstBranchAndBoundSolver, self).__init__()

    def search(self, initial):
        """
        (BreadthFirstBranchAndBoundSolver, Solution) -> Solution
        Does a breadth-first traversal of the solution space
        starting from the given node.
        """
        queue = QueueAsLinkedList()
        if initial.isFeasible:
            queue.enqueue(initial)
        while not queue.isEmpty:
            current = queue.dequeue()
            if current.isComplete:
                self.updateBest(current)
            else:
                for successor in current.successors:
                    if successor.isFeasible and \
                            successor.bound < self._bestObjective:
                        queue.enqueue(successor)
