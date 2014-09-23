#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: matrix.py,v $
#   $Revision: 1.25 $
#
#   $Id: matrix.py,v 1.25 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Matrix class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.25 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

#{
class Matrix(object):
    """
    Base class from which all matrices are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, numberOfRows, numberOfColumns):
        """
        (Matrix, int, int) -> None
        Constructs this matrix.
        """
	assert numberOfRows >= 0
	assert numberOfColumns >= 0
        super(Matrix, self).__init__()
        self._numberOfRows = numberOfRows
        self._numberOfColumns = numberOfColumns

    def getNumberOfRows(self):
        return self._numberOfRows

    numberOfRows = property(
        fget = lambda self: self.getNumberOfRows())

    def getNumberOfColumns(self):
        return self._numberOfColumns

    numberOfColumns = property(
        fget = lambda self: self.getNumberOfColumns())
#}>a

    transpose = property(
        fget = lambda self: self.getTranspose())

    @staticmethod
    def test(mat):
        "Matrix test program."
        print Matrix.test.__doc__
        k = 0
        for i in xrange(mat.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                mat[i, j] = k
                k = k + 1
        for i in xrange(mat.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                print mat[i, j],
            print
        mat = mat + mat
        for i in xrange(mat.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                print mat[i, j],
            print

    @staticmethod
    def testTranspose(mat):
        "Matrix transpose test program."
        print Matrix.testTranspose.__doc__
        mat[0,0] = 31
        mat[0,2] = 41
        mat[0,3] = 59
        mat[1,1] = 26
        mat[2,3] = 53
        mat[2,4] = 58
        mat[4,2] = 97
        mat[5,1] = 93
        mat[5,5] = 23
        for i in xrange(mat.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                print mat[i, j],
            print
        mat[2,4] = 0
        mat[5,3] = 0
        mat = mat.transpose
        for i in xrange(mat.numberOfRows):
            for j in xrange(mat.numberOfColumns):
                print mat[i, j],
            print

    @staticmethod
    def testTimes(mat1, mat2):
        "Matrix multiply test program."
        print Matrix.testTimes.__doc__
        mat1[0, 0] = 1
        mat1[0, 1] = 2
        mat1[0, 2] = 3
        mat2[0, 0] = 1
        mat2[1, 0] = 2
        mat2[2, 0] = 3
        for i in xrange(mat1.numberOfRows):
            for j in xrange(mat1.numberOfColumns):
                print mat1[i, j],
            print
        for i in xrange(mat2.numberOfRows):
            for j in xrange(mat2.numberOfColumns):
                print mat2[i, j],
            print
        mat1 = mat2 * mat1
        for i in xrange(mat1.numberOfRows):
            for j in xrange(mat1.numberOfColumns):
                print mat1[i, j],
            print
