#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: rectangle.py,v $
#   $Revision: 1.6 $
#
#   $Id: rectangle.py,v 1.6 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Rectangle class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.graphicalObject import GraphicalObject
from opus7.point import Point

#{
class Rectangle(GraphicalObject):
    """
    A rectangle.
    """

    def __init__(self, center, height, width):
        """
        (Rectangle, Point, int, int) -> None
        Constructs a rectangle with the given center, height, and width.
        """
        super(Rectangle, self).__init__(center)
        self._height = height
        self._width = width

    def draw(self):
        """
        (Rectangle) -> None
        Draws this rectangle.
        """
        # ...
#[
	print "RECTANGLE", self._center, self._height, self._width
#]
#}>a

    def _compareTo(self, c):
	raise MethodNotImplemented

    @staticmethod
    def main(*argv):
	"Rectangle test program."
	print Rectangle.main.__doc__
	c = Rectangle(Point(0,0), 1, 2)
	GraphicalObject.test(c)
	return 0

if __name__ == "__main__":
    sys.exit(Rectangle.main(*sys.argv))
