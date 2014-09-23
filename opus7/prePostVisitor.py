#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: prePostVisitor.py,v $
#   $Revision: 1.12 $
#
#   $Id: prePostVisitor.py,v 1.12 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PrePostVisitor class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.visitor import Visitor

#{
class PrePostVisitor(Visitor):
    """
    Pre/Post visitor class.
    """

#}@head

#{
#}@tail

#{
    def __init__(self):
        """
        (PrePostVisitor) -> None
        Constructor.
        """
        super(PrePostVisitor, self).__init__()

    def preVisit(self, obj):
        """
        (PrePostVisitor, Object) -> None
        Default pre-visit method does nothing.
        """
        pass

    def inVisit(self, obj):
        """
        (PrePostVisitor, Object) -> None
        Default in-visit method does nothing.
        """
        pass

    def postVisit(self, obj):
        """
        (PrePostVisitor, Object) -> None
        Default post-visit method does nothing.
        """
        pass

    visit = inVisit
#}>a

    @staticmethod
    def main(*argv):
        "PrePostVisitor test program."
        print PrePostVisitor.main.__doc__
        v = PrePostVisitor()
        print v.isDone
        return 0

if __name__ == "__main__":
    sys.exit(PrePostVisitor.main(*sys.argv))
