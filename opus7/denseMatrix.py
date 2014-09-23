#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: denseMatrix.py,v $
#   $Revision: 1.25 $
#
#   $Id: denseMatrix.py,v 1.25 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the DenseMatrix class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.25 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.matrix import Matrix
from opus7.multiDimensionalArray import MultiDimensionalArray

#{
class DenseMatrix(Matrix):
    """
    Dense matrix.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, rows, cols):
        """
        (DenseMatrix, int, int) -> None
        Constructor.
        """
        super(DenseMatrix, self).__init__(rows, cols)
        self._array = MultiDimensionalArray(rows, cols)

    def __getitem__(self, (i, j)):
        """
        (DenseMatrix, (int, int)) -> Object
        Returns the specified element of this matrix.
        """
        return self._array[i,j]

    def __setitem__(self, (i, j), value):
        """
        (DenseMatrix, (int, int), Object) -> Object
        Sets the specified element of this matrix to the given value.
        """
        self._array[i,j] = value
#}>a

    def getTranspose(self):
        """
        (DenseMatrix) -> DenseMatrix
        Returns the transpose of this matrix.
        """
        result = DenseMatrix(
            self.numberOfColumns, self.numberOfRows)
        for i in xrange(self.numberOfRows):
            for j in xrange(self.numberOfColumns):
                result[i,j] = self[j,i]
        return result

#{
    def __mul__(self, mat):
        """
        (DenseMatrix, DenseMatrix) -> DenseMatrix
        Returns the product of this matrix and the given matrix.
        """
        assert self.numberOfColumns == mat.numberOfRows
        result = DenseMatrix(
            self.numberOfRows, mat.numberOfColumns)
        for i in xrange(self.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                sum = 0
                for k in xrange(self.numberOfColumns):
                    sum += self[i,k] * mat[k,j]
                result[i,j] = sum
        return result
#}>b

    def __add__(self, mat):
        """
        (DenseMatrix, DenseMatrix) -> DenseMatrix
        Returns the sum of this matrix and the given matrix.
        """
        assert self.numberOfColumns == mat.numberOfColumns and \
            self.numberOfRows == mat.numberOfRows
        result = DenseMatrix(self.numberOfRows, self.numberOfColumns)
        for i in xrange(self.numberOfRows):
            for j in xrange(self.numberOfColumns):
                result[i,j] = self[i,j] + mat[i,j]
        return result

    @staticmethod
    def main(*argv):
        "DenseMatrix test program"
        print DenseMatrix.main.__doc__
        mat = DenseMatrix(6, 6)
        Matrix.test(mat)
        Matrix.testTranspose(mat)
        Matrix.testTimes(mat, mat)
        return 0

if __name__ == "__main__":
    sys.exit(DenseMatrix.main(*sys.argv))
