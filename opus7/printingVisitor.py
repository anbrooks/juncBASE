#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: printingVisitor.py,v $
#   $Revision: 1.6 $
#
#   $Id: printingVisitor.py,v 1.6 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the PrintingVisitor class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.visitor import Visitor

class PrintingVisitor(Visitor):
    """
    Visitor that prints the objects it visits.
    """

    def __init__(self):
        """
        (PrintingVisitor) -> None
        Constructor.
        """
        super(PrintingVisitor, self).__init__()
        self._comma = False

    def visit(self, obj):
        """
        (PrintingVisitor, Object) -> None
        Prints the object.
        """
        if self._comma:
            print ", ",
        print obj,
        self._comma = True

    def finish(self):
        """
        (PrintingVisitor) -> None
        Finishes the line.
        """
        print ""
        self._comma = False
