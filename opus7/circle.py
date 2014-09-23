#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: circle.py,v $
#   $Revision: 1.6 $
#
#   $Id: circle.py,v 1.6 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Circle class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.graphicalObject import GraphicalObject
from opus7.point import Point

#{
class Circle(GraphicalObject):
    """
    A circle.
    """

    def __init__(self, center, radius):
        """
        (Circle, Point, int) -> None
        Constructs a circle with the given center and radius.
        """
        super(Circle, self).__init__(center)
        self._radius = radius

    def draw(self):
        """
        (Circle) -> None
        Draws this circle.
        """
        # ...
#[
	print "CIRCLE", self._center,  self._radius
        pass
#]
#}>a

    def _compareTo(self, c):
	raise MethodNotImplemented

    @staticmethod
    def main(*argv):
	"Circle test program."
	print Circle.main.__doc__
	c = Circle(Point(0,0), 1)
	GraphicalObject.test(c)
	return 0

if __name__ == "__main__":
    sys.exit(Circle.main(*sys.argv))
