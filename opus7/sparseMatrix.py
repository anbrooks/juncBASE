#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sparseMatrix.py,v $
#   $Revision: 1.5 $
#
#   $Id: sparseMatrix.py,v 1.5 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SparseMatrix class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.5 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.matrix import Matrix

class SparseMatrix(Matrix):
    """
    Base class from which all sparse matrix classes are derived.
    """

    def __init__(self, numberOfRows, numberOfColums):
        """
        (SparseMatrix, int, int) -> None
        Constructs a sparse matrix with the given number of rows and columns.
        """
        super(SparseMatrix, self).__init__(numberOfRows, numberOfColums)
