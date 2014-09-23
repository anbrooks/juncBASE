#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: postOrder.py,v $
#   $Revision: 1.16 $
#
#   $Id: postOrder.py,v 1.16 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PostOrder class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.prePostVisitor import PrePostVisitor

#{
class PostOrder(PrePostVisitor):
    """
    Adapter to convert a Visitor to a PrePostVisitor for post-order traversal.
    """

    def __init__(self, visitor):
        """
        (PostOrder, Visitor) -> None
        Constructs a post-order visitor from the given visitor.
        """
        super(PostOrder, self).__init__()
        self._visitor = visitor

    def postVisit(self, obj):
        """
        (PostOrder, Object) -> None
        Post-visits the given object.
        """
        self._visitor.visit(obj)
    
    def getIsDone(self):
        """
        (PostOrder) -> bool
        Returns true if the visitor is done.
        """
        return self._visitor.isDone
#}>a

    @staticmethod
    def main(*argv):
        "PostOrder test program."
        print PostOrder.main.__doc__
        return 0

if __name__ == "__main__":
    sys.exit(PostOrder.main(*sys.argv))
