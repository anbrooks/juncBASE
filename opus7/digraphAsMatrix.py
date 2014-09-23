#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: digraphAsMatrix.py,v $
#   $Revision: 1.12 $
#
#   $Id: digraphAsMatrix.py,v 1.12 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the DigraphAsMatrix class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.iterator import Iterator
from opus7.digraph import Digraph
from opus7.graphAsMatrix import GraphAsMatrix

class DigraphAsMatrix(Digraph, GraphAsMatrix):
    """
    Directed graph implemented using and adjacency matrix.
    """

    def __init__(self, size):
        """
        (DigraphAsMatrix, int) -> None
        Constructs a digraph with the specified maximum number of vertices.
        """
        super(DigraphAsMatrix, self).__init__(size)

    def addEdge(self, *args):
        """
        (DigraphAsMatrxi, int, int [, Object]) -> None
        Adds an edge with the (optional) weight
        connecting the given vertices in this digraph.
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
        self._matrix [v, w] = self.Edge(self, v, w, weight)
        self._numberOfEdges += 1

    class EdgeIterator(Iterator):
        """
        Enumerates the edges of a DigraphAsMatrix.
        """

        def __init__(self, graph):
            """
            (DigraphAsMatrix.EdgeIterator, DigraphAsMatrix) -> None
            Constructs an iterator for the given digraph.
            """
            super(DigraphAsMatrix.EdgeIterator, self).__init__(graph)
            self._v = -1
            self._w = self._container._numberOfVertices - 1
        
        def next(self):
            """
            (DigraphAsMatrix.EdgeIterator) -> Edge
            Returns the next edge in this digraph.
            """
            while True:
                self._w += 1
                if self._w == self._container._numberOfVertices:
                    self._v += 1
                    self._w = 0
                    if self._v == self._container._numberOfVertices:
                        self._v = -1
                        self._w = self._container._numberOfVertices - 1
                        raise StopIteration
                if self._container._matrix[self._v, self._w] is not None:
                    break
            return self._container._matrix[self._v, self._w]

    def getEdges(self):
        """
        (DigraphAsMatrix) -> DigraphAsMatrix.EdgeIterator
        Returns an iterator that enumerates the edges in this digraph.
        """
        return self.EdgeIterator(self)
    
    @staticmethod
    def main(*argv):
        "DigraphAsMatrix test program."
        print DigraphAsMatrix.main.__doc__
        dg = DigraphAsMatrix(32)
        Digraph.test(dg)
        dg.purge()
        Digraph.testWeighted(dg)
        return 0

if __name__ == "__main__":
    sys.exit(DigraphAsMatrix.main(*sys.argv))
