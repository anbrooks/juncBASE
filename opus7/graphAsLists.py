#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: graphAsLists.py,v $
#   $Revision: 1.18 $
#
#   $Id: graphAsLists.py,v 1.18 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the GraphAsLists class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.18 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.graph import Graph
from opus7.array import Array
from opus7.linkedList import LinkedList
from opus7.iterator import Iterator

#{
class GraphAsLists(Graph):
    "Graph implemented using adjacency lists."

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size):
        """
        (GraphAsLists, int) -> None
        Constructs a graph with the given maximum number of vertices.
        """
        super(GraphAsLists, self).__init__(size)
        self._adjacencyList = Array(size)
        for i in xrange(size):
            self._adjacencyList[i] = LinkedList()
#}>a

    def purge(self):
        """
        (GraphAsLists) -> None
        Purges this graph.
        """
        for i in xrange(self._numberOfVertices):
            self._adjacencyList[i].purge()
        super(GraphAsLists, self).purge()

    def addEdge(self, *args):
        """
        (GraphAsLists, int, int [, Object]) -> None
        Adds an edge with the (optional) given weight
        connecting the given vertices in this graph.
        """
        if len(args) == 2:
            v = args[0]
            w = args[1]
            weight = None
        elif len(args) == 3:
            v = args[0]
            w = args[1]
            weight = args[2]
        else:
            raise ValueError
        self._adjacencyList[v].append(self.Edge(self, v, w, weight))
        self._numberOfEdges += 1

    def getEdge(self, v, w):
        """
        (GraphAsLists, int, int) -> Edge
        Returns the edge connecting the specified vertices in this graph.
        """
        if v < 0 or v >= self._numberOfVertices:
            raise IndexError
        if w < 0 or w >= self._numberOfVertices:
            raise IndexError
        ptr = self._adjacencyList[v].head
        while ptr is not None:
            edge = ptr.datum
            if edge.v1.number == w:
                return edge
            ptr = ptr.next
        raise KeyError

    def isEdge(self, v, w):
        """
        (GraphAsLists, int, int) -> bool
        Returns true if there is an edge connecting the specified vertices
        in this graph.
        """
        if v < 0 or v >= self._numberOfVertices:
            raise IndexError
        if w < 0 or w >= self._numberOfVertices:
            raise IndexError
        ptr = self._adjacencyList[v].head
        while ptr is not None:
            edge = ptr.datum
            if edge.v1.number == w:
                return True
            ptr = ptr.next
        return False

    class EdgeIterator(Iterator):
        """
        Enumerates the edges of a GraphAsLists.
        """

        def __init__(self, graph):
            """
            (GraphAsLists.EdgeIterator, GraphAsLists) -> None
            Constructs an iterator for the given graph.
            """
            super(GraphAsLists.EdgeIterator, self).__init__(graph)
            self._v = -1
            self._ptr = None

        def next(self):
            """
            (GraphAsLists.EdgeIterator) -> Edge
            Returns the next edge in this graph.
            """
            if self._ptr is not None:
                self._ptr = self._ptr.next
            if self._ptr is None:
                self._v += 1
                while self._v < self._container._numberOfVertices:
                    self._ptr = self._container._adjacencyList[self._v].head
                    if self._ptr is not None:
                        break
                    self._v += 1
            if self._ptr is None:
                self._v = -1
                raise StopIteration
            return self._ptr.datum

    def getEdges(self):
        """
        (GraphAsLists) -> GraphAsLists.EdgeIterator
        Returns an iterator that enumerates the edges in this graph.
        """
        return self.EdgeIterator(self)

    class EmanatingEdgeIterator(Iterator):
        """
        Enumerates the edges emanating from a given vertex in this graph.
        """

        def __init__(self, graph, v):
            """
            (GraphAsLists.EmanatingEdgesIterator, GraphAsLists, int) -> None
            Constructs an iterator that enumerates the edges
            emanating from the given vertex of the given graph.
            """
            super(GraphAsLists.EmanatingEdgeIterator, self).__init__(graph)
            self._v = v
            self._ptr = None

        def next(self):
            """
            (GraphAsLists.EmanatingEdgesIterator) -> Edge
            Returns the next edge emanating from the vertex in the graph.
            """
            if self._ptr is None:
                self._ptr = self._container._adjacencyList[self._v].head
            else:
                self._ptr = self._ptr.next
            if self._ptr is None:
                raise StopIteration
            return self._ptr.datum

    def getEmanatingEdges(self, v):
        """
        (GraphAsLists, int) -> GraphAsLists.EmanatingEdgesIterator
        Returns an iterator that enumerates the edges emanating from
        the given vertex in this graph.
        """
        return self.EmanatingEdgeIterator(self, v)

    def getIncidentEdges(self, v):
        """
        (GraphAsLists, int) -> GraphAsLists.EmanatingEdgesIterator
        Returns an iterator that enumerates the edges incident upon
        the given vertex in this graph.
        """
        raise NotImplementedError

    def getIsCyclic(self):
        """
        (GraphAsLists) -> bool

        Returns true if this graph is cyclic.
        """
        raise NotImplementedError

    def _compareTo(self, obj):
        """
        (GraphAsLists, GraphAsLists) -> int

        Compares this graph with the given graph.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "GraphAsLists test program."
        print GraphAsLists.main.__doc__
        g = GraphAsLists(32)
        Graph.test(g)
        g.purge()
        Graph.testWeighted(g)
        return 0

if __name__ == "__main__":
    sys.exit(GraphAsLists.main(*sys.argv))
