#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sparseMatrixAsLinkedList.py,v $
#   $Revision: 1.12 $
#
#   $Id: sparseMatrixAsLinkedList.py,v 1.12 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SparseMatrixAsLinkedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.matrix import Matrix
from opus7.sparseMatrix import SparseMatrix
from opus7.array import Array
from opus7.linkedList import LinkedList

class SparseMatrixAsLinkedList(SparseMatrix):
    """
    Sparse matrix implemented using a linked list of matrix entries.
    """

    def __init__(self, numberOfRows, numberOfColumns):
        """
        (SparseMatrixAsLinkedList, int, int) -> None
        Constructs a sparse matrix with the given number of rows and columns.
        """
        super(SparseMatrixAsLinkedList, self).__init__(
            numberOfRows, numberOfColumns)
        self._lists = Array(numberOfRows)
        for i in xrange(numberOfRows):
            self._lists[i] = LinkedList()

    def getNumberOfRows(self):
        """
        (SparseMatrixAsLinkedList) -> int
        Returns the number of rows in this matrix.
        """
        return self._numberOfRows

    def getNumberOfColumns(self):
        """
        (SparseMatrixAsLinkedList) -> int
        Returns the number of columns in this matrix.
        """
        return self._numberOfColumns

    def __getitem__(self, indices):
        """
        (SparseMatrixAsLinkedList, (int, int)) -> Object
        Returns the element of this matrix at the given indices.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        ptr = self._lists[i].head
        while ptr is not None:
            entry = ptr.datum
            if entry._column == j:
                return entry._datum
            if entry._column > j:
                break
            ptr = ptr.next
        return 0

    def __setitem__(self, indices, value):
        """
        (SparseMatrixAsLinkedList, (int, int), Object) -> None
        Sets the element of this matrix at the given indices to the given value.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        ptr = self._lists[i].head
        while ptr is not None:
            entry = ptr.datum
            if entry._column == j:
                entry._datum = value
                return
            elif entry._column > j:
                ptr.insertBefore(self.Entry(i, j, value))
                return
            ptr = ptr.next
        self._lists[i].append(self.Entry(i, j, value))

    def putZero(self, i, j):
        """
        (SparseMatrixAsLinkedList, (int, int), Object) -> None
        Sets the element of this matrix at the given indices to zero.
        """
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        ptr = lists[i].head
        while ptr is not None:
            entry = ptr.datum
            if entry._column == j:
                self._lists[i].extract(entry)
                return
            ptr = ptr.next

    def getTranspose(self):
        """
        (SparseMatrixAsLinkedList) -> SparseMatrixAsLinkedList
        Returns the transpose of this matrix.
        """
        result = SparseMatrixAsLinkedList(
            self._numberOfColumns, self._numberOfRows)
        for i in xrange(self._numberOfColumns):
            ptr = self._lists[i].head
            while ptr is not None:
                entry = ptr.datum
                result._lists[entry._column].append(
                    self.Entry(entry._column, entry._row, entry._datum))
                ptr = ptr.next
        return result

    class Entry(object):
        """
        Represents a element of a sparse matrix.
        """

        def __init__(self, row, column, datum):
            """
            (SparseMatrixAsLinkedList.Entry, int, int, Object) -> None
            Constructs an element of the sparse matrix
            with the given indices and value.
            """
            self._row = row
            self._column = column
            self._datum = datum

    @staticmethod
    def main(*argv):
        "SparseMatrixAsLinkedList test program."
        print SparseMatrixAsLinkedList.main.__doc__
        mat = SparseMatrixAsLinkedList(6, 6)
        #Matrix.TestMatrix(mat)
        Matrix.testTranspose(mat)
        #Matrix.TestTimes(mat, mat)
        return 0

if __name__ == "__main__":
    sys.exit(SparseMatrixAsLinkedList.main(*sys.argv))
