#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/12/07 01:32:17 $
#   $RCSfile: algorithms.py,v $
#   $Revision: 1.26 $
#
#   $Id: algorithms.py,v 1.26 2005/12/07 01:32:17 brpreiss Exp $
#

"""
Provides the Algorithms class.
"""

__author__ = "Bruno R. Preiss, P.Eng."
__date__ = "$Date: 2005/12/07 01:32:17 $"
__version__ = "$Revision: 1.26 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object
from opus7.array import Array
from opus7.denseMatrix import DenseMatrix
from opus7.stackAsLinkedList import StackAsLinkedList
from opus7.queueAsLinkedList import QueueAsLinkedList
from opus7.chainedHashTable import ChainedHashTable
from opus7.association import Association
from opus7.avlTree import AVLTree
from opus7.partitionAsForest import PartitionAsForest
from opus7.binaryHeap import BinaryHeap
from opus7.graphAsLists import GraphAsLists
from opus7.digraphAsMatrix import DigraphAsMatrix
from opus7.digraphAsLists import DigraphAsLists
from opus7.postOrder import PostOrder
from opus7.visitor import Visitor

#{
class Algorithms(object):
    """
    Contains a bunch of algorithms.
    """

#}@head

#{
#}@tail

#{
    @staticmethod
    def breadthFirstTraversal(tree):
        """
        (Tree) -> None
        Does a breadth-first traversal of a tree and prints the keys.
        """
        queue = QueueAsLinkedList()
#!        queue.enqueue(tree)
#[
	if not tree.isEmpty:
	    queue.enqueue(tree)
#]
        while not queue.isEmpty:
            t = queue.dequeue()
            print t.key
            for i in xrange(t.degree):
                subTree = t.getSubtree(i)
#!		queue.enqueue(subTree)
#[
		if not subTree.isEmpty:
		    queue.enqueue(subTree)
#]
#}>a

#{
    @staticmethod
    def equivalenceClasses(input, output):
        """
        (File, File) -> None
        Computes equivalence classes using a partition.
        First reads an integer from the input stream that
        specifies the size of the universal set.
        Then reads pairs of integers from the input stream
        that denote equivalent items in the universal set.
        Prints the partition on end-of-file.
        """
        line = input.readline()
        p = PartitionAsForest(int(line))
        for line in input.readlines():
            words = line.split()
            i = int(words[0])
            j = int(words[1])
            s = p.find(i)
            t = p.find(j)
            if s is not t:
                p.join(s, t)
            else:
                output.write("redundant pair: %d, %d\n" % (i, j))
        output.write(str(p) + "\n")
#}>b

#{
    @staticmethod
    def DijkstrasAlgorithm(g, s):
        """
        (Digraph, int) -> DigraphAsLists
        Dijkstra's algorithm to solve the single-source, shortest path problem
        for the given edge-weighted, directed graph.
        """
        n = g.numberOfVertices
        table = Array(n)
        for v in xrange(n):
            table[v] = Algorithms.Entry()
        table[s].distance = 0
        queue = BinaryHeap(g.numberOfEdges)
        queue.enqueue(Association(0, g[s]))
        while not queue.isEmpty:
            assoc = queue.dequeueMin()
            v0 = assoc.value
            if not table[v0.number].known:
                table[v0.number].known = True
                for e in v0.emanatingEdges:
                    v1 = e.mateOf(v0)
                    d = table[v0.number].distance + e.weight
                    if table[v1.number].distance > d:

                        table[v1.number].distance = d
                        table[v1.number].predecessor = v0.number
                        queue.enqueue(Association(d, v1))
        result = DigraphAsLists(n)
        for v in xrange(n):
            result.addVertex(v, table[v].distance)
        for v in xrange(n):
            if v != s:
                result.addEdge(v, table[v].predecessor)
        return result
#}>c

#{
    @staticmethod
    def FloydsAlgorithm(g):
        """
        (Digraph) -> DigraphAsMatrix
        Floyd's algorithm to solve the all-pairs, shortest path problem
        for the given edge-weighted, directed graph.
        """
        n = g.numberOfVertices
        distance = DenseMatrix(n, n)
        for v in xrange(n):
            for w in xrange(n):
                distance[v, w] = sys.maxint
        for e in g.edges:
            distance[e.v0.number, e.v1.number] = e.weight
        for i in xrange(n):
            for v in xrange(n):
                for w in xrange(n):
                    if distance[v, i] != sys.maxint and \
                        distance[i, w] != sys.maxint:
                        d = distance[v, i] + distance[i, w]
                        if distance[v, w] > d:
                            distance[v, w] = d
        result = DigraphAsMatrix(n)
        for v in xrange(n):
            result.addVertex(v)
        for v in xrange(n):
            for w in xrange(n):
                if distance[v, w] != sys.maxint:
                    result.addEdge(v, w, distance[v, w])
        return result
#}>d

#{
    @staticmethod
    def PrimsAlgorithm(g, s):
        """
        (Graph, int) -> GraphAsLists
        Prim's algorithm to find a minimum-cost spanning tree
        for the given edge-weighted, undirected graph.
        """
        n = g.numberOfVertices
        table = Array(n)
        for v in xrange(n):
            table[v] = Algorithms.Entry()
        table[s].distance = 0
        queue = BinaryHeap(g.numberOfEdges)
        queue.enqueue(Association(0, g[s]))
        while not queue.isEmpty:
            assoc = queue.dequeueMin()
            v0 = assoc.value
            if not table[v0.number].known:
                table[v0.number].known = True
                for e in v0.emanatingEdges:
                    v1 = e.mateOf(v0)
                    d = e.weight
                    if not table[v1.number].known and \
                            table[v1.number].distance > d:
                        table[v1.number].distance = d
                        table[v1.number].predecessor = v0.number
                        queue.enqueue(Association(d, v1))
        result = GraphAsLists(n)
        for v in xrange(n):
            result.addVertex(v)
        for v in xrange(n):
            if v != s:
                result.addEdge(v, table[v].predecessor)
        return result
#}>e

#{
    @staticmethod
    def KruskalsAlgorithm(g):
        """
        (Graph) -> GraphAsLists
        Kruskal's algorithm to find a minimum-cost spanning tree
        for the given edge-weighted, undirected graph.
        """
        n = g.numberOfVertices
        result = GraphAsLists(n)
        for v in xrange(n):
            result.addVertex(v)
        queue = BinaryHeap(g.numberOfEdges)
        for e in g.edges:
            weight = e.weight
            queue.enqueue(Association(weight, e))
        partition = PartitionAsForest(n)
        while not queue.isEmpty and partition.count > 1:
            assoc = queue.dequeueMin()
            e = assoc.value
            n0 = e.v0.number
            n1 = e.v1.number
            s = partition.find(n0)
            t = partition.find(n1)
            if s != t:
                partition.join(s, t)
                result.addEdge(n0, n1)
        return result
#}>f

#{
    class Entry(object):
        """
        Data structure used in Dijkstra's and Prim's algorithms.
        """

        def __init__(self):
            """
            (Algorithms.Entry) -> None
            Constructor.
            """
            self.known = False
            self.distance = sys.maxint
            self.predecessor = sys.maxint
#}>g

#{
    @staticmethod
    def criticalPathAnalysis(g):
        """
        (Digraph) -> DigraphAsLists
        Computes the critical path in an event-node graph.
        """
        n = g.numberOfVertices

        earliestTime = Array(n)
        earliestTime[0] = 0
        g.topologicalOrderTraversal(
            Algorithms.EarliestTimeVisitor(earliestTime))

        latestTime = Array(n)
        latestTime[n - 1] = earliestTime[n - 1]
        g.depthFirstTraversal(PostOrder(
            Algorithms.LatestTimeVisitor(latestTime)), 0)

        slackGraph = DigraphAsLists(n)
        for v in xrange(n):
            slackGraph.addVertex(v)
        for e in g.edges:
            slack = latestTime[e.v1.number] - \
                earliestTime[e.v0.number] - e.weight
            slackGraph.addEdge(
                e.v0.number, e.v1.number, slack)
        return Algorithms.DijkstrasAlgorithm(slackGraph, 0)
#}>h

#{
    class EarliestTimeVisitor(Visitor):
        """
        Used by the critical path analysis program
        to compute the earliest completion time for each event.
        """

        def __init__(self, earliestTime):
            """
            (Algorithms.EarliestTimeVisitor, Array) -> None
            Constructor.
            """
            super(Algorithms.EarliestTimeVisitor,self).__init__()
            self._earliestTime = earliestTime

        def visit(self, w):
            """
            (Algorithms.EarliestTimeVisitor, Vertex) -> None
            Determines the earliest completion time for the given vertex.
            """
            t = self._earliestTime[0]
            for e in w.incidentEdges:
                t = max(t,
                    self._earliestTime[e.v0.number] + e.weight)
            self._earliestTime[w.number] = t
#}>i

    class LatestTimeVisitor(Visitor):
        """
        Used by the critical path analysis program
        to compute the latest completion time for each event.
        """

        def __init__(self, latestTime):
            """
            (Algorithms.LatestTimeVisitor, Array) -> None
            Constructor.
            """
            super(Algorithms.LatestTimeVisitor, self).__init__()
            self._latestTime = latestTime

        def visit(self, v):
            """
            (Algorithms.LatestTimeVisitor, Vertex) -> None
            Determines the latest completion time for the given vertex.
            """
            t = self._latestTime[len(self._latestTime) - 1]
            for e in v.emanatingEdges:
                t = min(t,
                    self._latestTime[e.v1.number] - e.weight)
            self._latestTime[v.number] = t

#{
    @staticmethod
    def calculator(input, output):
        """
        (File, File) -> None
        A very simple reverse-Polish calculator.
        """
        stack = StackAsLinkedList()
        for line in input.readlines():
            for word in line.split():
                if word == "+":
                    arg2 = stack.pop()
                    arg1 = stack.pop()
                    stack.push(arg1 + arg2)
                elif word == "*":
                    arg2 = stack.pop()
                    arg1 = stack.pop()
                    stack.push (arg1 * arg2)
                elif word == "=":
                    arg = stack.pop()
                    output.write(str(arg) + "\n")
                else:
                    stack.push(int(word))
#}>j

#{
    class Counter(object):
        """
        A counter.
        """

        def __init__(self, value):
            """
            (Algorithms.Counter, int) -> None
            Constructs a counter with the given initial value.
            """
            super(Algorithms.Counter, self).__init__()
            self._value = value

        def __repr__(self):
            """
            (Algorithms.Counter) -> str
            Returns a string representation of this counter.
            """
            return str(self._value)

        def __iadd__(self, value):
            """
            (Algorithms.Counter, int) -> None
	    Add the give value to this counter.
            """
            self._value += value

    @staticmethod
    def wordCounter(input, output):
        """
        (File, File) -> None
        Counts the number of occurrences of each word in the given file.
        """
        table = ChainedHashTable(1031)
        for line in input.readlines():
            for word in line.split():
                assoc = table.find(Association(word))
                if assoc is None:
                    table.insert(Association(
                        word, Algorithms.Counter(1)))
                else:
                    counter = assoc.value
		    counter += 1
        output.write(str(table) + "\n")
#}>k

#{
    @staticmethod
    def translate(dictionary, input, output):
        """
        (File, File, File) -> None
        Reads all the word pairs from the dictionary file
        and then reads words from the input file,
        translates the words (if possible),
        and writes them to the output file.
        """
        searchTree = AVLTree()
        for line in dictionary.readlines():
            words = line.split()
	    assert len(words) == 2
            searchTree.insert(Association(words[0], words[1]))
        for line in input.readlines():
            for word in line.split():
                assoc = searchTree.find(Association(word))
                if assoc is None:
                    output.write(word + " ")
                else:
                    output.write(assoc.value + " ")
            output.write("\n")
#}>l
