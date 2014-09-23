#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: digraph.py,v $
#   $Revision: 1.21 $
#
#   $Id: digraph.py,v 1.21 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Digraph class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.21 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.graph import Graph
from opus7.array import Array
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.printingVisitor import PrintingVisitor
from opus7.preOrder import PreOrder

#{
class Digraph(Graph):
    """
    Base class from which all directed graph classes are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size):
        """
        (Digraph, int) -> None
        Constructs a digraph with the specified maximum number of vertices.
        """
        super(Digraph, self).__init__(size)
        self._isDirected = True

    @abstractmethod
    def getIsStronglyConnected(self): pass

    isStronglyConnected = property(
        fget = lambda self: self.getIsStronglyConnected())

    @abstractmethod
    def topologicalOrderTraversal(self): pass
#}>a

#{
    def topologicalOrderTraversal(self, visitor):
        """
        (Graph, Visitor) -> None
        Makes the given visitor visit the vertices of this graph
        in topological order.
        """
        inDegree = Array(self.numberOfVertices)
        for v in xrange(self.numberOfVertices):
            inDegree[v] = 0
        for e in self.edges:
            inDegree[e.v1.number] += 1
        queue = QueueAsLinkedList()
        for v in xrange(self.numberOfVertices):
            if inDegree[v] == 0:
                queue.enqueue(self[v])
        while not queue.isEmpty and not visitor.isDone:
            v = queue.dequeue()
            visitor.visit(v)
            for to in v.successors:
                inDegree[to.number] -= 1
                if inDegree[to.number] == 0:
                    queue.enqueue(to)
#}>b

#{
    def getIsStronglyConnected(self):
        """
        (Graph) -> bool
        Returns true if the graph is strongly connected.
        """
        for v in xrange(self.numberOfVertices):
            visitor = self.CountingVisitor()
            self.depthFirstTraversal(PreOrder(visitor), v)
            if visitor.count != self.numberOfVertices:
                return False
        return True
#}>c

#{
    def getIsCyclic(self):
        """
        (Graph) -> bool
        Returns true if the graph is cyclic.
        """
        visitor = self.CountingVisitor()
        self.topologicalOrderTraversal(visitor)
        return visitor.count != self.numberOfVertices
#}>d

    @staticmethod
    def test(g):
        "Digraph test program."
        print Digraph.test.__doc__
        Graph.test(g)
        visitor = PrintingVisitor()
        print "TopologicalOrderTraversal"
        g.topologicalOrderTraversal(visitor)
        visitor.finish()
        print "isCyclic returns " + str(g.isCyclic)
        print "isStronglyConnected returns " + str(g.isStronglyConnected)

    @staticmethod
    def testWeighted(g):
        "Weighted digraph test program."
        print Digraph.testWeighted.__doc__
        Graph.testWeighted(g)
