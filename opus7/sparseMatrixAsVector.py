#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: sparseMatrixAsVector.py,v $
#   $Revision: 1.14 $
#
#   $Id: sparseMatrixAsVector.py,v 1.14 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the SparseMatrixAsVector class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.14 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.matrix import Matrix
from opus7.sparseMatrix import SparseMatrix
from opus7.array import Array

class SparseMatrixAsVector(SparseMatrix):
    """
    Sparse matrix implemented as a vector of non-zero entries.
    """

    def __init__(self, numberOfRows, numberOfColumns, numberOfElements):
        """
        (SparseMatrixAsVector, int, int, int) -> None
        Constructs a sparse matrix with the given number of rows and columns
        and the given number of non-zero entries.
        """
        super(SparseMatrixAsVector, self).__init__(
            numberOfRows, numberOfColumns)
        self._numberOfElements = numberOfElements
        self._array = Array(numberOfElements)
        for i in xrange(numberOfElements):
            self._array[i] = self.Entry()

    def getNumberOfRows(self):
        """
        (SparseMatrixAsVector) -> int
        Returns the number of rows in this matrix.
        """
        return self._numberOfRows

    def getNumberOfColumns(self):
        """
        (SparseMatrixAsVector) -> int
        Returns the number of columns in this matrix.
        """
        return self._numberOfColumns

    def findPosition(self, i, j):
        """
        (SparseMatrixAsVector, int, int) -> int
        Returns the position in the vector of the entry with the given indices.
        """
        target = i * self._numberOfColumns + j
        left = 0
        right = self._numberOfElements - 1
        while left <= right:
            middle = (left + right) / 2
            probe = self._array[middle]._row * self._numberOfColumns \
                + self._array[middle]._column
            if target > probe:
                left = middle + 1
            elif target < probe:
                right = middle - 1
            else:
                return middle
        return -1

    def __getitem__(self, indices):
        """
        (SparseMatrixAsVector, (int, int)) -> Object
        Returns the element of this matrix at the given indices.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = self.findPosition(i, j)
        if position >= 0:
            return self._array[position]._datum
        else:
            return 0

    def __setitem__(self, indices, value):
        """
        (SparseMatrixAsVector, (int, int), Object) -> None
        Sets the element of this matrix at the given indices to the given value.
        """
        i = indices[0]
        j = indices[1]
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = self.findPosition(i, j)
        if position >= 0:
            self._array[position]._datum = value
        else:
            if len(self._array) == self._numberOfElements:
                newArray = Array(2 * len(self._array))
                for p in xrange(len(self._array)):
                    newArray[p] = self._array[p]
                for p in xrange(len(self._array), len(newArray)):
                    newArray[p] = self.Entry()
                self._array = newArray
            k = self._numberOfElements
            while k > 0 and (self._array[k - 1]._row > i or \
                    self._array[k - 1]._row == i and \
                    self._array[k - 1]._column >= j):
                self._array[k] = self._array[k - 1]
                k -= 1
            self._array[k] = self.Entry(i, j, value)
            self._numberOfElements += 1

    def putZero(self, i, j):
        """
        (SparseMatrixAsVector, int, int) -> None
        Sets the element of this matrix at the given indices to zero.
        """
        if i < 0 or i >= self._numberOfRows:
            raise IndexError
        if j < 0 or j >= self._numberOfColumns:
            raise IndexError
        position = self.findPosition(i, j)
        if position >= 0:
            self._numberOfElements -= 1
            for k in xrange(position, self._numberOfElements):
                self._array[k] = self._array[k + 1]
            self._array[k] = self.Entry()

    def getTranspose(self):
        """
        (SparseMatrixAsVector) -> SparseMatrixAsVector
        Returns the transpose of this matrix.
        """
        result  = SparseMatrixAsVector(
            self._numberOfColumns, self._numberOfRows, self._numberOfElements)
        offset = Array(self._numberOfColumns)
        for i in xrange(self._numberOfColumns):
            offset[i] = 0
        for i in xrange(self._numberOfElements):
            offset[self._array[i]._column] += 1
        sum = 0
        for i in xrange(self._numberOfColumns):
            tmp = offset[i]
            offset[i] = sum
            sum += tmp
        for i in xrange(self._numberOfElements):
            result._array[offset[self._array[i]._column]] = \
                self.Entry(self._array[i]._column, self._array[i]._row,
                    self._array[i]._datum)
            offset[self._array[i]._column] += 1
        result._numberOfElements = self._numberOfElements
        return result

    def __mul__(self, mat):
        """
        (SparseMatrixAsVector, SparseMatrixAsVector) -> SparseMatrixAsVector
        Returns the product of this matrix and the given matrix.
        """
        assert self._numberOfColumns == mat._numberOfRows
        matT = mat.transpose
        result = SparseMatrixAsVector(
            self._numberOfRows, matT.numberOfRows,
            self._numberOfRows + matT.numberOfRows)
        iPosition = 0
        while iPosition < self._numberOfElements:
            i = self._array[iPosition]._row
            jPosition = 0
            while jPosition < matT._numberOfElements:
                j = matT._array[jPosition]._row
                sum = 0
                k1 = iPosition
                k2 = jPosition
                while k1 < self._numberOfElements \
                        and self._array[k1]._row == i \
                        and k2 < matT._numberOfElements \
                        and matT._array[k2]._row == j:
                    if self._array[k1]._column < matT._array[k2]._column:
                        k1 += 1
                    elif self._array[k1]._column > matT._array[k2]._column:
                        k2 += 1
                    else:
                        sum += self._array[k1]._datum * matT._array[k2]._datum
                        k1 += 1
                        k2 += 1
                if sum != 0:
                    result[i, j] = sum
                while jPosition < matT._numberOfElements and \
                        matT._array[jPosition]._row == j:
                    jPosition += 1
            while iPosition < self._numberOfElements and \
                    self._array[iPosition]._row == i:
                iPosition += 1
        return result

    class Entry(object):
        """
        Represents an entry in this sparse matrix.
        """

        def __init__(self, row=0, column=0, datum=0):
            """
            (SparseMatrixAsVector.Entry, int, int, Object) -> None
            Constructs an entry with the given indices and value.
            """
            self._row = row
            self._column = column
            self._datum = datum

    @staticmethod
    def main(*argv):
        "SparseMatrixAsVector test program."
        print SparseMatrixAsVector.main.__doc__
        mat = SparseMatrixAsVector(6, 6, 12)
        #Matrix.TestMatrix(mat)
        Matrix.testTranspose(mat)
        Matrix.testTimes(mat, mat)
        return 0

if __name__ == "__main__":
    sys.exit(SparseMatrixAsVector.main(*sys.argv))
