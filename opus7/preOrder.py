#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: preOrder.py,v $
#   $Revision: 1.16 $
#
#   $Id: preOrder.py,v 1.16 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PreOrder class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.prePostVisitor import PrePostVisitor

#{
class PreOrder(PrePostVisitor):
    """
    Adapter to convert a Visitor to a PrePostVisitor for pre-order traversal.
    """

    def __init__(self, visitor):
        """
        (PreOrder, Visitor) -> None
        Constructs a pre-order visitor from the given visitor.
        """
        super(PreOrder, self).__init__()
        self._visitor = visitor

    def preVisit(self, obj):
        """
        (PreOrder, Object) -> None
        Pre-visits the given object.
        """
        self._visitor.visit(obj)
    
    def getIsDone(self):
        """
        (PreOrder) -> bool
        Returns true if the visitor is done.
        """
        return self._visitor.isDone
#}>a

    @staticmethod
    def main(*argv):
        "PreOrder test program."
        print PreOrder.main.__doc__
        return 0

if __name__ == "__main__":
    sys.exit(PreOrder.main(*sys.argv))
