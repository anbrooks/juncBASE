#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:37 $
#   $RCSfile: association.py,v $
#   $Revision: 1.25 $
#
#   $Id: association.py,v 1.25 2005/06/09 00:00:37 brpreiss Exp $
#

"""
Provides the Association class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:37 $"
__version__ = "$Revision: 1.25 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object

#{
class Association(Object):
    """
    Represents a (key, value) pair using a tuple.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, *args):
        """
        (Association, Object [,Object]) -> None
        Constructs an association with the given key and optional value.
        """
        if len(args) == 1:
            self._tuple = (args[0], None)
        elif len(args) == 2:
            self._tuple = args
        else:
            raise ValueError
#}>a

#{
    def getKey(self):
        """
        (Association) -> Object
        Returns the key of this association.
        """
        return self._tuple[0]

    key = property(
        fget = lambda self: self.getKey())

    def getValue(self):
        """
        (Association) -> Object
        Returns the value of this association.
        """
        return self._tuple[1]

    value = property(
        fget = lambda self: self.getValue())
#}>b

#{
    def _compareTo(self, assoc):
        """
        (Association, Assocation) -> int
        Compares the key of this association
        with with the key of given association.
        """
        assert isinstance(self, assoc.__class__)
        return cmp(self.key, assoc.key)

    def __str__(self):
        """
        (Association) -> string
        Returns a string representation of this association.
        """
        return "Association %s" % str(self._tuple)
#}>c

#{
    def __hash__(self):
        """
        (Association) -> int
        Hashes the key of this association.
        """
        return hash(self.key)
#}>d

    @staticmethod
    def main(*argv):
        "Association test program. "
        print Association.main.__doc__
        a = Association(1,2)
        print a, a.key, a.value
        print hash(a)
        print Association(2,2) > Association(3,2)
        b = Association(3)
        print b
        return 0

if __name__ == "__main__":
    sys.exit(Association.main(*sys.argv))
