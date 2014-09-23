#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sparseMatrixAsArray.py,v $
#   $Revision: 1.12 $
#
#   $Id: sparseMatrixAsArray.py,v 1.12 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SparseMatrixAsArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.matrix import Matrix
from opus7.sparseMatrix import SparseMatrix
from opus7.multiDimensionalArray import MultiDimensionalArray

class SparseMatrixAsArray(SparseMatrix):
    """
    Sparse matrix implemented using a two-dimensional array.
    """

    END_OF_ROW = -1

    def __init__(self, numberOfRows, numberOfColumns, fill):
        """
        (SparseMatrixAsArray, int, int, int) -> None
        Constructs a sparse matrix with the given number of rows and columns
        and with the given row fill (maximum number of items in a row).
        """
        super(SparseMatrixAsArray, self).__init__(numberOfRows, numberOfColumns)
        self._fill = fill
        self._values = MultiDimensionalArray(numberOfColumns, fill)
        self._columns = MultiDimensionalArray(numberOfColumns, fill)
        for i in xrange(self._numberOfRows):
            self._columns[i, 0] = self.END_OF_ROW

    def getNumberOfRows(self):
        """
        (SparseMatrixAsArray) -> int
        Returns the number of rows in this matrix.
        """
        return self._numberOfRows

    def getNumberOfColumns(self):
        """
        (SparseMatrixAsArray) -> int
        Returns the number of columns in this matrix.
        """
        return self._numberOfColumns

    def findPosition(self, i, j):
        """
        (SparseMatrixAsArray, int, int) -> int
        Returns the column in which the (i,j) element
        of this matrix is stored.
        """
        k = 0
        while k < self._fill and self._columns[i, k] != self.END_OF_ROW:
            if self._columns[i, k] == j:
                return k
            k += 1
        return -1

    def __getitem__(self, indices):
        """
        (SparseMatrixAsArray, (int, int)) -> Object
        Returns the object at the given indices in this array.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = self.findPosition(i, j)
        if position >= 0:
            return self._values[i, position]
        else:
            return 0

    def __setitem__(self, indices, value):
        """
        (SparseMatrixAsArray, (int, int)) -> Object
        Sets the object at the given indices in this array to the given value.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            print i, j
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = self.findPosition(i, j)
        if position >= 0:
            self._values[i, position] = value
        else:
            k = 0
            while k < self._fill and self._columns[i, k] != self.END_OF_ROW:
                k += 1
            if k >= self._fill:
                raise IndexError("row is full")
            if k < self._fill - 1:
                self._columns[i, k + 1] = self.END_OF_ROW
            while k > 0 and self._columns[i, k] >= j:
                self._values[i, k] = self._values[i, k - 1]
                self._columns[i, k] = self._columns[i, k - 1]
                k -= 1
            self._values[i, k] = value
            self._columns[i, k] = j

    def putZero(self, i, j):
        """
        (SparseMatrixAsArray, (int, int)) -> None
        Zeroes the element at the given indices in this array.
        """
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = findPosition(i, j)
        if position >= 0:
            k = position
            while k < numberOfColumns - 1 and \
                    self._columns[i, k + 1] != self.END_OF_ROW:
                self._values[i, k] = self._values[i, k + 1]
                self._columns[i, k] = self._columns[i, k + 1]
                k += 1
            if k < self._numberOfColumns:
                self._columns[i, k] = self.END_OF_ROW

    @staticmethod
    def main(*argv):
        "SparseMatrixAsArray test program."
        print SparseMatrixAsArray.main.__doc__
        mat = SparseMatrixAsArray(6, 6, 3)
        #Matrix.test(mat)
        #Matrix.TestTranspose(mat)
        #Matrix.TestTimes(mat, mat)
        return 0

if __name__ == "__main__":
    sys.exit(SparseMatrixAsArray.main(*sys.argv))
