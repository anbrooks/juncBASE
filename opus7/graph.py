#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: graph.py,v $
#   $Revision: 1.29 $
#
#   $Id: graph.py,v 1.29 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Graph class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.29 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import types
from opus7.abstractmethod import abstractmethod
from opus7.container import Container
from opus7.vertex import Vertex
from opus7.edge import Edge
from opus7.iterator import Iterator
from opus7.visitor import Visitor
from opus7.prePostVisitor import PrePostVisitor
from opus7.printingVisitor import PrintingVisitor
from opus7.preOrder import PreOrder
from opus7.array import Array
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.exception import *

#{
class Graph(Container):
    """
    Base class from which all graph classes are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, size):
        """
        (Graph, int) -> None
        Constructs a graph with the given maximum number of vertices.
        """
        super(Graph, self).__init__()
        self._numberOfVertices = 0
        self._numberOfEdges = 0
        self._vertex = Array(size)
        self._isDirected = False

    class Vertex(Vertex):
        """
        Represents a vertex in a graph.
        """

        def __init__(self, graph, number, weight):
            """
            (Graph.Vertex, Graph, int, Object) -> None
            Constructs a vertex of the given graph
            with the given number and weight.
            """
            super(Graph.Vertex, self).__init__()
            self._graph = graph
            self._number = number
            self._weight = weight

#!        # ...
#[
        def getNumber(self):
            """
            (Graph.Vertex) -> int
            Returns the number of this vertex.
            """
            return self._number

        def getWeight(self):
            """
            (Graph.Vertex) -> Object
            Returns the number of this vertex.
            """
            return self._weight
        
        def __hash__(self):
            """
            (Graph.Vertex) -> int
            Returns the hash of this vertex.
            """
            result = self._number
            if self._weight is not None:
                result += hash(self._weight)
            return result
        
        def __str__(self):
            """
            (Graph.Vertex) -> str
            Returns a string representation of this vertex.
            """
            if self._weight is None:
                return "Vertex {%d}" % (self._number)
            else:
                return "Vertex {%d, weight = %s}" % (
                    self._number, str(self._weight))

        def getIncidentEdges(self):
            """
            (Graph.Vertex) -> iterator
            Returns an iterator that enumerates
            the incident edges of this vertex.
            """
            return self._graph.getIncidentEdges(self._number)

        def getEmanatingEdges(self):
            """
            (Graph.Vertex) -> iterator
            Returns an iterator that enumerates
            the emanating edges of this vertex.
            """
            return self._graph.getEmanatingEdges(self._number)
        
        class PredecessorIterator(Iterator):
            """
            Enumerates the predecessors of this vertex.
            """

            def __init__(self, vertex):
                """
                (Graph.Vertex.PredecessorIterator, Graph.Vertex) -> None
                Constructs an iterator that enumerates
                the predecessors of the given vertex.
                """
                super(Graph.Vertex.PredecessorIterator, self).__init__(
                    self._vertex._graph)
                self._vertex = vertex
                self._edges = vertex.incidentEdges

            def next(self):
                """
                (Graph.Vertex.PredecessorIterator) -> Vertex
                Returns the next predecessor.
                """
                return self._edges.next().mateOf(self._vertex)
            
        def getPredecessors(self):
            """
            (Graph.Vertex) -> Graph.Vertex.PredecessorIterator
            Returns an iterator that enumerates the predecesors of this vertex.
            """
            return self.PredecessorIterator(self)

        class SuccessorIterator(Iterator):
            """
            Enumerates the successors of this vertex.
            """

            def __init__(self, vertex):
                """
                (Graph.Vertex.SuccessorIterator, Graph.Vertex) -> None
                Constructs an iterator that enumerates
                the successors of the given vertex.
                """
                super(Graph.Vertex.SuccessorIterator, self).__init__(
                    vertex._graph)
                self._vertex = vertex
                self._edges = vertex.emanatingEdges

            def next(self):
                """
                (Graph.Vertex.SuccessorIterator) -> Vertex
                Returns the next successor.
                """
                return self._edges.next().mateOf(self._vertex)
            
        def getSuccessors(self):
            """
            (Graph.Vertex) -> Graph.Vertex.SuccessorIterator
            Returns an iterator that enumerates the successors of this vertex.
            """
            return self.SuccessorIterator(self)

        def _compareTo(self, obj):
            """
            (Graph.Vertex, Graph.Vertex) -> int

            Compares this vertex with the given vertex.
            """
            assert isinstance(self, obj.__class__)
            raise NotImplementedError
#]

    class Edge(Edge):
        """
        Represents an edge in a graph.
        """

        def __init__(self, graph, v0, v1, weight):
            """
            (Graph.Edge, Graph, int, int, Object) -> None
            Constructs an edge of the given graph
            between the given vertices and with the given weight.
            """
            super(Graph.Edge, self).__init__()
            self._graph = graph
            self._v0 = v0
            self._v1 = v1
            self._weight = weight

#!        # ...
#[
        def getV0(self):
            """
            (Graph.Edge) -> Vertex
            Returns the first vertex of this edge.
            """
            return self._graph._vertex[self._v0]

        def getV1(self):
            """
            (Graph.Edge) -> Vertex
            Returns the second vertex of this edge.
            """
            return self._graph._vertex[self._v1]
        
        def getWeight(self):
            """
            (Graph.Edge) -> Weight
            Returns the weight of this edge.
            """
            return self._weight
        
        def mateOf(self, v):
            """
            (Graph.Edge, Vertex) -> Vertex
            Returns the mate of the given vertex of this edge.
            """
            if v.number == self._v0:
                return self._graph._vertex[self._v1]
            elif v.number == self._v1:
                return self._graph._vertex[self._v0]
            else:
                raise ValueError

        def getIsDirected(self):
            """
            (Graph.Edge) -> bool
            Returns true if this edge is directed.
            """
            return self._graph.isDirected
        
        def __hash__(self):
            """
            (Graph.Edge) -> int
            Returns the hash of this edge.
            """
            result = self._v0 * self._graph._numberOfVertices + self._v1
            if self._weight is not None:
                result += hash(self._weight)
            return result
        
        def __str__(self):
            """
            (Graph.Edge) -> int
            Returns the string representation of this edge.
            """
            s = "Edge {" + str(self._v0)
            if self.isDirected:
                s = s + "->" + str(self._v1)
            else:
                s = s + "--" + str(self._v1)
            if self._weight is not None:
                s = s + ", weight = " + str(self._weight)
            s = s + "}"
            return s

        def _compareTo(self, obj):
            """
            (Graph.Edge, Graph.Edge) -> int

            Compares this edge with the given edge.
            """
            raise NotImplementedError
#]
#}>a

#{
    @abstractmethod
    def getNumberOfEdges(self): pass

    numberOfEdges = property(
        fget = lambda self: self.getNumberOfEdges())

    @abstractmethod
    def getNumberOfVertices(self): pass

    numberOfVertices = property(
        fget = lambda self: self.getNumberOfVertices())

    @abstractmethod
    def getIsDirected(self): pass

    isDirected = property(
        fget = lambda self: self.getIsDirected())

    @abstractmethod
    def getIsConnected(self): pass

    isConnected = property(
        fget = lambda self: self.getIsConnected())

    @abstractmethod
    def getIsCyclic(self): pass

    isCyclic = property(
        fget = lambda self: self.getIsCyclic())

    @abstractmethod
    def getVertices(self): pass

    vertices = property(
        fget = lambda self: self.getVertices())

    @abstractmethod
    def getEdges(self): pass

    edges = property(
        fget = lambda self: self.getEdges())
#}>b

#{
    @abstractmethod
    def addVertex(self, *args): pass

    @abstractmethod
    def getVertex(self, v): pass

    @abstractmethod
    def addEdge(self, *args): pass

    @abstractmethod
    def getEdge(self, v, w): pass

    @abstractmethod
    def isEdge(self, v, w): pass

    @abstractmethod
    def depthFirstTraversal(self, visitor, start): pass

    @abstractmethod
    def breadthFirstTraversal(self, visitor, start): pass

    @abstractmethod
    def getIncidentEdges(self, v): pass

    @abstractmethod
    def getEmanatingEdges(self, v): pass

    def __len__(self):
        return self.numberOfVertices

    def __getitem__(self, v):
        return self.getVertex(v)
#}>c

    def purge(self):
        """
        (Graph) -> None
        Purges this graph.
        """
        for i in xrange(self._numberOfVertices):
            self._vertex[i] = None
        self._numberOfVertices = 0
        self._numberOfEdges = 0

    def addVertex(self, *args):
        """
        (Graph, int [, Object]) -> None
        Adds a vertex (with optional weight) to this graph.
        """
        if self._numberOfVertices == len(self._vertex):
            raise ContainerFull
        if len(args) == 1:
            v = args[0]
            weight = None
        elif len(args) == 2:
            v = args[0]
            weight = args[1]
        else:
            raise ValueError
        if args[0] != self._numberOfVertices:
            raise ValueError
        self._vertex[self._numberOfVertices] = self.Vertex(self, v, weight)
        self._numberOfVertices += 1

    def getVertex(self, v):
        """
        (Graph, int) -> Vertex
        Returns the specified vertex of this graph.
        """
        if v < 0 or v >= self._numberOfVertices:
            raise IndexError
        return self._vertex[v]

    def getIsDirected(self):
        """
        (Graph) -> bool
        Returns true if this graph is a directed graph.
        """
        return self._isDirected

    def getNumberOfVertices(self):
        """
        (Graph) -> int
        Returns the number of vertices in this graph.
        """
        return self._numberOfVertices

    def getNumberOfEdges(self):
        """
        (Graph) -> int
        Returns the number of edges in this graph.
        """
        return self._numberOfEdges

    def accept(self, visitor):
        """
        (Graph, Visitor) -> None
        Makes the given visitor visit all the vertices in this graph.
        """
        assert isinstance(visitor, Visitor)
        for v in xrange(self._numberOfVertices):
            if visitor.isDone:
                break
            visitor.visit(self._vertex[v])

    class Iterator(Iterator):
        """
        Enumerates the vertices of a graph.
        """

        def __init__(self, graph):
            """
            (Graph.Iterator, Graph) -> None
            Constructs an iterator for the given graph.
            """
            super(Graph.Iterator, self).__init__(graph)
            self._v = -1

        def next(self):
            """
            (Graph.Iterator) -> Vertex
            Returns the next vertex in this graph.
            """
            self._v += 1
            if self._v == self._container._numberOfVertices:
                self._v = -1
                raise StopIteration
            return self._container._vertex[self._v]

    def getVertices(self):
        """
        (Graph) -> Graph.Iterator
        Returns an iterator for this graph.
        """
        return self.Iterator(self)

#{
    def depthFirstTraversal(self, visitor, start):
        """
        (Graph, PrePostVisitor, Vertex) -> None
        Makes the given visitor visit the vertices of this graph
        in depth-first traversal order starting from the given vertex.
        """
        assert isinstance(visitor, PrePostVisitor)
        visited = Array(self._numberOfVertices)
        for v in xrange(self._numberOfVertices):
            visited[v] = False
        self._depthFirstTraversal(visitor, self[start], visited)

    def _depthFirstTraversal(self, visitor, v, visited):
        """
        (Graph, PrePostVisitor, Vertex, Array) -> None
        Implements the depth-first traversal
        using an array to keep track of already visited vertices.
        """
        if visitor.isDone:
            return
        visitor.preVisit(v)
        visited[v.number] = True
        for to in v.successors:
            if not visited[to.number]:
                self._depthFirstTraversal(visitor, to, visited)
        visitor.postVisit(v)
#}>d

#{
    def breadthFirstTraversal(self, visitor, start):
        """
        (Graph, Visitor, Vertex) -> None
        Makes the given visitor visit the vertices of this graph
        in breadth-first traversal order starting from the given vertex.
        """
        assert isinstance(visitor, Visitor)
        enqueued = Array(self._numberOfVertices)
        for v in xrange(self._numberOfVertices):
            enqueued[v] = False
        queue = QueueAsLinkedList()
        queue.enqueue(self[start])
        enqueued[start] = True
        while not queue.isEmpty and not visitor.isDone:
            v = queue.dequeue()
            visitor.visit(v)
            for to in v.successors:
                if not enqueued[to.number]:
                    queue.enqueue(to)
                    enqueued[to.number] = True
#}>e

    class StrVisitor(Container.StrVisitor):
        """
        Visitor that accumuates visited items in a string.
        """

        def __init__(self):
            """
            (Graph.StrVisitor) -> None
            Constructor.
            """
            super(Graph.StrVisitor, self).__init__()

        def visit(self, v):
            """
            (Graph.StrVisitor, Vertex) -> None
            Appends the given vertex and all the edges emanating
            from that vertex to the string.
            """
            self._string = self._string + str(v) + "\n"
            for e in v.emanatingEdges:
                self._string = self._string + "    " + str(e) + "\n"

    def __str__(self):
        """
        (Graph) -> str
        Returns the string representation of this graph.
        """
        visitor = self.StrVisitor()
        self.accept(visitor)
        return self.__class__.__name__ + " {\n" + str(visitor) + "}"

#{
    class CountingVisitor(Visitor):
        """
        A visitor that counts the objects it visits.
        """

        def __init__(self):
            """
            (Graph.CountingVisitor) -> None
            Constructs a counting visitor.
            """
            super(Graph.CountingVisitor, self).__init__()
            self._count = 0

        def visit(self, obj):
            """
            (Graph.CountingVisitor, Object) -> None
            Visits the given object.
            """
            self._count += 1
        
        def getCount(self):
            return self._count

        count = property(
            fget = lambda self: self.getCount())

    def getIsConnected(self):
        """
        (Graph) -> bool
        Returns true if the graph is connected.
        """
        visitor = self.CountingVisitor()
        self.depthFirstTraversal(PreOrder(visitor), 0)
        return visitor.count == self.numberOfVertices
#}>f

    def __iter__(self):
        """
        (Graph) -> Graph.Iterator
        Returns an iterator that enumerates the vertices of this graph.
        """
        return self.vertices

    def test(g):
        "Graph test program."
        print Graph.test.__doc__
        g.addVertex(0)
        g.addVertex(1)
        g.addVertex(2)
        g.addEdge(0, 1)
        g.addEdge(0, 2)
        g.addEdge(1, 2)
        print g
        print "isDirected returns " + str(g.isDirected)
        print "Using vertex iterator"
        for v in g.vertices:
            print v
        print "Using edge iterator"
        for e in g.edges:
            print e
        visitor = PrintingVisitor()
        print "DepthFirstTraversal"
        g.depthFirstTraversal(PreOrder(visitor), 0)
        visitor.finish()
        print "BreadthFirstTraversal"
        g.breadthFirstTraversal(visitor, 0)
        visitor.finish()
        print "isConnected returns " + str(g.isConnected)

    def testWeighted(g):
        "Weighted graph test program."
        print Graph.testWeighted.__doc__
        g.addVertex(0, 123)
        g.addVertex(1, 234)
        g.addVertex(2, 345)
        g.addEdge(0, 1, 3)
        g.addEdge(0, 2, 1)
        g.addEdge(1, 2, 4)
        print g
        print "Using vertex iterator"
        for v in g.vertices:
            print v
        print "Using edge iterator"
        for e in g.edges:
            print e
