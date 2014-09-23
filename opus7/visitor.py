#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: visitor.py,v $
#   $Revision: 1.26 $
#
#   $Id: visitor.py,v 1.26 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the Visitor class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.26 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.object import Object

#{
class Visitor(Object):
    """
    Visitor class.
    """

#}@head

#{
#}@tail

#{
    def __init__(self):
        """
        (Visitor) -> None
        Constructs this visitor.
        """
        super(Visitor, self).__init__()

    def visit(self, obj):
        """
        (Visitor, Object) -> None
        Default visit method does nothing.
        """
        pass
    
    def getIsDone(self):
        """
        (Visitor) -> bool
        Default isDone_get method returns false always.
        """
        return False

    isDone = property(
        fget = lambda self: self.getIsDone())
#}>a

    def _compareTo(self, obj):
        """
        (Visitor, Visitor) -> int

        Compares this visitor with the given visitor.
        """
        assert isinstance(self, obj.__class__)
        raise NotImplementedError

    @staticmethod
    def main(*argv):
        "Visitor test program."
        print Visitor.main.__doc__
        v = Visitor()
        print v.isDone
        return 0

if __name__ == "__main__":
    sys.exit(Visitor.main(*sys.argv))
