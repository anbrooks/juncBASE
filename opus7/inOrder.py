#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: inOrder.py,v $
#   $Revision: 1.16 $
#
#   $Id: inOrder.py,v 1.16 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the InOrder class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.prePostVisitor import PrePostVisitor

#{
class InOrder(PrePostVisitor):
    """
    Adapter to convert a Visitor to a PrePostVisitor for in-order traversal.
    """

    def __init__(self, visitor):
        """
        (InOrder, Visitor) -> None
        Constructs a in-order visitor from the given visitor.
        """
        super(InOrder, self).__init__()
        self._visitor = visitor

    def inVisit(self, obj):
        """
        (InOrder, Object) -> None
        In-visits the given object.
        """
        self._visitor.visit(obj)
    
    def getIsDone(self):
        """
        (InOrder) -> bool
        Returns true if the visitor is done.
        """
        return self._visitor.isDone
#}>a

    @staticmethod
    def main(*argv):
        "InOrder test program."
        print InOrder.main.__doc__
        return 0

if __name__ == "__main__":
    sys.exit(InOrder.main(*sys.argv))
