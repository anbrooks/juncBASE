#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: graphAsMatrix.py,v $
#   $Revision: 1.18 $
#
#   $Id: graphAsMatrix.py,v 1.18 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the GraphAsMatrix class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.18 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.graph import Graph
from opus7.denseMatrix import DenseMatrix
from opus7.iterator import Iterator

#{
class GraphAsMatrix(Graph):
    """
    Graph implemented using an adjacency matrix.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size):
        """
        (GraphAsMatrix, int) -> None
        Constructs a graph with the given maximum number of vertices.
        """
        super(GraphAsMatrix, self).__init__(size)
        self._matrix = DenseMatrix(size, size)
#}>a

    def purge(self):
        """
        (GraphAsMatrix) -> None
        Purges this graph.
        """
        for i in xrange(self._numberOfVertices):
            for j in xrange(self._numberOfVertices):
                self._matrix[i, j] = None
        super(GraphAsMatrix, self).purge()

    def addEdge(self, *args):
        """
        (GraphAsMatrix, int, int [, Object]) -> None
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
        if self._matrix[v, w] is not None:
            raise ValueError
        if v == w:
            raise ValueError
        edge = self.Edge(self, v, w, weight)
        self._matrix[v, w] = edge
        self._matrix[w, v] = edge
        self._numberOfEdges += 1

    def getEdge(self, v, w):
        """
        (GraphAsMatrix, int, int) -> Edge
        Returns the edge connecting the specified vertices in this graph.
        """
        if v < 0 or v >= self._numberOfVertices:
            raise IndexError
        if w < 0 or w >= self._numberOfVertices:
            raise IndexError
        if self._matrix[v, w] is None:
            raise KeyError
        return self._matrix[v, w]

    def isEdge(self, v, w):
        """
        (GraphAsMatrix, int, int) -> bool
        Returns true if there is an edge connecting the specified vertices
        in this graph.
        """
        if v < 0 or v >= self._numberOfVertices:
            raise IndexError
        if w < 0 or w >= self._numberOfVertices:
            raise IndexError
        return matrix[v, w] is not None

    class EdgeIterator(Iterator):
        """
        Enumerates the edges of a GraphAsMatrix.
        """

        def __init__(self, graph):
            """
            (GraphAsMatrix.EdgeIterator, GraphAsMatrix) -> None
            Constructs an iterator for the given graph.
            """
            super(GraphAsMatrix.EdgeIterator, self).__init__(graph)
            self._v = -1
            self._w = self._container._numberOfVertices - 1
        
        def next(self):
            """
            (GraphAsMatrix.EdgeIterator) -> Edge
            Returns the next edge in this graph.
            """
            while True:
                self._w += 1
                if self._w == self._container._numberOfVertices:
                    self._v += 1
                    self._w = self._v
                    if self._v == self._container._numberOfVertices:
                        self._v = -1
                        self._w = self._container._numberOfVertices - 1
                        raise StopIteration
                if self._container._matrix[self._v, self._w] is not None:
                    break
            return self._container._matrix[self._v, self._w]

    def getEdges(self):
        """
        (GraphAsMatrix) -> GraphAsMatrix.EdgeIterator
        Returns an iterator that enumerates the edges in this graph.
        """
        return self.EdgeIterator(self)

    class EmanatingEdgesIterator(Iterator):
        """
        Enumerates the edges emanating from a given vertex in this graph.
        """

        def __init__(self, graph, v):
            """
            (GraphAsMatrix.EmanatingEdgesIterator, GraphAsMatrix, int) -> None
            Constructs an iterator that enumerates the edges
            emanating from the given vertex of the given graph.
            """
            super(GraphAsMatrix.EmanatingEdgesIterator, self).__init__(graph)
            self._v = v
            self._w = -1

        def next(self):
            """
            (GraphAsMatrix.EmanatingEdgesIterator) -> Edge
            Returns the next edge emanating from the vertex in the graph.
            """
            while True:
                self._w += 1
                if self._w == self._container._numberOfVertices:
                    self._w = -1
                    raise StopIteration
                if self._container._matrix[self._v, self._w] is not None:
                    break
            return self._container._matrix[self._v, self._w]

    def getEmanatingEdges(self, v):
        """
        (GraphAsMatrix, int) -> GraphAsMatrix.EmanatingEdgesIterator
        Returns an iterator that enumerates the edges emanating from
        the given vertex in this graph.
        """
        return self.EmanatingEdgesIterator(self, v)

    class IncidentEdgesIterator(Iterator):
        """
        Enumerates the edges incident upon a given vertex in this graph.
        """

        def __init__(self, graph, w):
            """
            (GraphAsMatrix.IncidentEdgesIterator, GraphAsMatrix, int) -> None
            Constructs an iterator that enumerates the edges
            incident upon the given vertex of the given graph.
            """
            super(GraphAsMatrix.IncidentEdgesIterator, self).__init__(graph)
            self._v = -1
            self._w = w

        def next(self):
            """
            (GraphAsMatrix.IncidentEdgesIterator) -> Edge
            Returns the next edge incident upon the vertex in the graph.
            """
            while True:
                self._v += 1
                if self._v == self._container._numberOfVertices:
                    self._v = -1
                    raise StopIteration
                if self._container._matrix[self._v, self._w] is not None:
                    break
            return self._container._matrix[self._v, self._w]

    def getIncidentEdges(self, w):
        """
        (GraphAsMatrix, int) -> GraphAsMatrix.IncidentEdgesIterator
        Returns an iterator that enumerates the edges incident upon
        the given vertex in this graph.
        """
        return self.IncidentEdgesIterator(self, w)

    def getIsCyclic(self):
        """
        (GraphAsMatrix) -> bool

        Returns true if this graph is cyclic.
        """
        raise NotImplementedError

    def _compareTo(self, obj):
        """
        (GraphAsMatrix, GraphAsMatrix) -> int

        Compares this graph with the given graph.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "GraphAsMatrix test program."
        print GraphAsMatrix.main.__doc__
        g = GraphAsMatrix(32)
        Graph.test(g)
        g.purge()
        Graph.testWeighted(g)
        return 0

if __name__ == "__main__":
    sys.exit(GraphAsMatrix.main(*sys.argv))
