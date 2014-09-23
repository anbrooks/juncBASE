#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: square.py,v $
#   $Revision: 1.6 $
#
#   $Id: square.py,v 1.6 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the Square class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.rectangle import Rectangle
from opus7.graphicalObject import GraphicalObject
from opus7.point import Point

#{
class Square(Rectangle):
    """
    A square.
    """

    def __init__(self, center, width):
        """
        (Square, Point, width) -> None
        Constructs a square with the given width (and height).
        """
        super(Square, self).__init__(center, width, width)
#}>a

    def _compareTo(self, c):
	raise MethodNotImplemented

    @staticmethod
    def main(*argv):
	"Square test program."
	print Square.main.__doc__
	c = Square(Point(0,0), 1)
	GraphicalObject.test(c)
	return 0

if __name__ == "__main__":
    sys.exit(Square.main(*sys.argv))
