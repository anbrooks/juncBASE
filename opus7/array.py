#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:37 $
#   $RCSfile: array.py,v $
#   $Revision: 1.33 $
#
#   $Id: array.py,v 1.33 2005/06/09 00:00:37 brpreiss Exp $
#

"""
Provides the Array class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:37 $"
__version__ = "$Revision: 1.33 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from exceptions import IndexError

#{
class Array(object):
    """
    Array class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, length = 0, baseIndex = 0):
        """
        (Array, int) -> None
        Constructs an array of the given length.
        """
	assert length >= 0
        self._data = [ None for i in xrange(length) ]
        self._baseIndex = baseIndex
#}>a

#{
    def __copy__(self):
        """
        (Array) -> Array
        Returns a shallow copy of this array.
        """
        result = Array(len(self._data))
        for i, datum in enumerate(self._data):
            result._data[i] = datum
        result._baseIndex = self._baseIndex
        return result
#}>b

#{
    def getOffset(self, index):
        """
        (Array, int) -> int
        Returns the offset for the given index.
        """
        offset = index - self._baseIndex
        if offset < 0 or offset >= len(self._data):
            raise IndexError
        return offset

    def __getitem__(self, index):
        """
        (Array, int) -> Object
        Returns the item in this array at the given index.
        """
        return self._data[self.getOffset(index)]

    def __setitem__(self, index, value):
        """
        (Array, int, Object) -> None
        Sets the item in this array at the given index to the given value.
        """
        self._data[self.getOffset(index)] = value
#}>c

#{
    def getData(self):
        return self._data

    data = property(
        fget = lambda self: self.getData())

    def getBaseIndex(self):
        return self._baseIndex

    def setBaseIndex(self, baseIndex):
        self._baseIndex = baseIndex

    baseIndex = property(
        fget = lambda self: self.getBaseIndex(),
        fset = lambda self, value: self.setBaseIndex(value))
#}>d

#{
    def __len__(self):
        """
        (Array) -> int
        Returns the length of this array.
        """
        return len(self._data)

    def setLength(self, value):
        if len(self._data) != value:
            newData = [ None for i in xrange(value) ]
            m = min(len(self._data), value)
            for i in xrange(m):
                newData[i] = self._data[i]
            self._data = newData

    length = property(
        fget = lambda self: self.__len__(),
        fset = lambda self, value: self.setLength(value))
#}>e

    def __str__(self):
        """
        (Array) -> string
        Returns a string representation of this array.
        """
        return "Array {baseIndex = %d, data = %s}" % (
            self._baseIndex, str(self._data))

    @staticmethod
    def main(*argv):
        "Array test program. "
        print Array.main.__doc__
        a1 = Array(3)
        print "a1 = %s" % (a1)
        a1[0] = 2
        a1[1] = a1[0] + 2
        a1[2] = a1[1] + 2
        print "a1 = %s" % (a1)
        print "baseIndex = ", a1.baseIndex
        print "length = ", a1.length
        a2 = Array(2, 10)
        a2[10] = 57
        print "a2 = %s" % (a2)
        print "baseIndex = ", a2.baseIndex
        print "length = ", a2.length
        return 0

if __name__ == "__main__":
    sys.exit(Array.main(*sys.argv))
