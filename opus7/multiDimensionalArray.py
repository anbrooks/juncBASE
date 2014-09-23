#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: multiDimensionalArray.py,v $
#   $Revision: 1.13 $
#
#   $Id: multiDimensionalArray.py,v 1.13 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the MultiDimensionalArray class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.13 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.array import Array

#{
class MultiDimensionalArray(object):
    """
    Multi-dimensional array implemented using a one-dimensional array.
    """

#}@head

#{

    # ...
#}@tail


#{
    def __init__(self, *dimensions):
        """
        (MultiDimensionalArray, Array) -> None
        Constructs a multi-dimensional array with the given dimensions.
        """
        self._dimensions = Array(len(dimensions))
        self._factors = Array(len(dimensions))
        product = 1
        i = len(dimensions) - 1
        while i >= 0:
            self._dimensions[i] = dimensions[i]
            self._factors[i] = product
            product *= self._dimensions[i]
            i -= 1
        self._data = Array(product)
#}>a

#{
    def getOffset(self, indices):
        """
        (MultiDimensionalArray, Array) -> int
        Maps the given indices of this multi-dimensional array
        into a one-dimensional array index.
        """
        if len(indices) != len(self._dimensions):
            raise IndexError
        offset = 0
        for i, dim in enumerate(self._dimensions):
            if indices[i] < 0 or indices[i] >= dim:
                raise IndexError
            offset += self._factors[i] * indices[i]
        return offset

    def __getitem__(self, indices):
        """
        (MultiDimensionalArray, Array) -> Object
        Returns the object in this multi-dimensional array
        at the given indices.
        """
        return self._data[self.getOffset(indices)]

    def __setitem__(self, indices, value):
        """
        (MultiDimensionalArray, Array) -> Object
        Sets the object in this multi-dimensional array
        at the given indices to the given value.
        """
        self._data[self.getOffset(indices)] = value
#}>b

    def __str__(self):
        """
        (MultiDimensionalArray) -> str
        Returns a string representation of this multi-dimensional array.
        """
        return \
            "MultiDimensionalArray {dimensions = %s, factors = %s, data = %s}" \
            % (str(self._dimensions), str(self._factors), str(self._data))

    @staticmethod
    def main(*argv):
        "MultiDimensionalArray test program."
        print MultiDimensionalArray.main.__doc__
        m = MultiDimensionalArray(2, 3, 4)
        m[1, 2, 3] = 57
        print m[1, 2, 3]
        print m
        return 0

if __name__ == "__main__":
    sys.exit(MultiDimensionalArray.main(*sys.argv))
