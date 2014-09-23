#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: integer.py,v $
#   $Revision: 1.6 $
#
#   $Id: integer.py,v 1.6 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Integer class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import warnings
from opus7.object import Object

#{
class Integer(int, Object):
    """
    Integer class.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, obj):
        """
        (Integer, object) -> None
        Constructs a string with the string representation of the given object.
        The Object metaclass provides a __new__ method
        that initializes the str instance, so none is defined here.
        """
        pass

    def _compareTo(self, obj):
        """
        (Integer, Integer) -> int

        Compares this string with the given string.
        """
        assert isinstance(self, obj.__class__)
        return cmp(int(self), int(obj))

#{
    def __hash__(self):
        """
        (Integer) -> int
        Hashes this string.
        """
	return self & sys.maxint
#}>a

    @staticmethod
    def testHash():
        "Integer hash test program."
        print Integer.testHash.__doc__
	print "57=%d" % hash(Integer(57))
	print "-123=%d" % hash(Integer(-123))

    @staticmethod
    def main(*argv):
        "Integer test program."
        print Integer.main.__doc__
        Integer.testHash()
        return 0

if __name__ == "__main__":
    sys.exit(Integer.main(*sys.argv))
