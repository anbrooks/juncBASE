#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sorter.py,v $
#   $Revision: 1.22 $
#
#   $Id: sorter.py,v 1.22 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Sorter class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.22 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.object import Object
from opus7.array import Array
from opus7.randomNumberGenerator import RandomNumberGenerator
from opus7.timer import Timer

#{
class Sorter(Object):
    """
    Base class from which all sorters are derived.
    """

    def __init__(self):
        """
        (Sorter) -> None
        Constructor.
        """
        super(Sorter, self).__init__()
        self._array = None
        self._n = 0

    @abstractmethod
    def _sort(self): pass

    def sort(self, array):
        """
        (Sorter, Array) -> None
        Sorts the given array.
        """
        assert isinstance(array, Array)
        self._array = array
        self._n = len(array)
        if self._n > 0:
            self._sort()
        self._array = None

    def swap(self, i, j):
        """
        (Sorter, int, int) -> None
        Swaps the specified element of the array being sorted.
        """
        tmp = self._array[i]
        self._array[i] = self._array[j]
        self._array[j] = tmp
#}>a

    def _compareTo(self, obj):
        """
        (Sorter, Sorter) -> int

        Compares this sorter with the given sorter.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def test(sorter, n, seed, *args):
        if len(args) == 0:
            m = 0
        elif len(args) == 1:
            m = args[0]
        else:
            raise ValueError
        RandomNumberGenerator.seed = seed
        data = Array(n)
        for i in xrange(n):
            datum = int(sys.maxint * RandomNumberGenerator.next)
            if m != 0:
                datum %= m
            data[i] = datum
        timer = Timer()
        timer.start()
        sorter.sort(data)
        timer.stop()
        datum = "%s %s %s %g" % (
            sorter.__class__.__name__, n, seed, timer.getElapsedTime())
        print datum
        sys.stderr.write(datum + "\n")
        for i in xrange(1, n):
            if data[i] < data[i - 1]:
                print "FAILED"
                break
