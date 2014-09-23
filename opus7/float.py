#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: float.py,v $
#   $Revision: 1.6 $
#
#   $Id: float.py,v 1.6 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Float class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import math
from opus7.object import Object

#{
class Float(float, Object):
    """
    Float class.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, obj):
        """
        (Float, object) -> None
        Constructs a string with the string representation of the given object.
        The Object metaclass provides a __new__ method
        that initializes the str instance, so none is defined here.
        """
        pass

    def _compareTo(self, obj):
        """
        (Float, Float) -> float

        Compares this string with the given string.
        """
        assert isinstance(self, obj.__class__)
        return cmp(float(self), float(obj))

#{
    def __hash__(self):
        """
        (Float) -> float
        Hashes this string.
        """
	(m, e) = math.frexp(self)
	mprime = int((abs(m) - 0.5) * (1L << 52))
	return mprime >> 20
#}>a

    @staticmethod
    def testHash():
        "Float hash test program."
        print Float.testHash.__doc__
	print "57.=0%o" % hash(Float(57.))
	print "23.=0%o" % hash(Float(23.))
	print "0.75=0%o" % hash(Float(0.75))
	print "-123.e6=0%o" % hash(Float(-123.e6))
	print "-123.e7=0%o" % hash(Float(-123.e7))
	print "0.875=0%o" % hash(Float(0.875))
	print "14.=0%o" % hash(Float(14.))

    @staticmethod
    def main(*argv):
        "Float test program."
        print Float.main.__doc__
        Float.testHash()
        return 0

if __name__ == "__main__":
    sys.exit(Float.main(*sys.argv))
