#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: graphicalObject.py,v $
#   $Revision: 1.7 $
#
#   $Id: graphicalObject.py,v 1.7 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the GraphicalObject class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.7 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.object import Object
from opus7.abstractmethod import abstractmethod
from opus7.point import Point

#{
class GraphicalObject(Object):
    """
    Base class from which all graphical objects are derived.
    """

#}@head

#{
#}@tail

    BACKGROUND_COLOR = 0
    FOREGROUND_COLOR = 1

    def setPenColor(self, color):
        pass

#{
    def __init__(self, center):
        """
        (GraphicalObject, Point) -> None
        Constructs a graphical object with the given point as its center.
        """
        super(GraphicalObject, self).__init__()
        self._center = center

    @abstractmethod
    def draw(self):
        """
        (GraphicalObject) -> None
        Draws this graphical object.
        """
        pass

    def erase(self):
        """
        (GraphicalObject) -> None
        Erases this graphical object.
        """
        self.setPenColor(self.BACKGROUND_COLOR)
        self.draw()
        self.setPenColor(self.FOREGROUND_COLOR)

    def moveTo(self, p):
        """
        (GraphicalObject) -> None
        Moves the center of this graphical object to the given point.
        """
        self.erase()
        self._center = p
        self.draw()
#}>a

    @staticmethod
    def test(go):
	"GraphicalObject test program."
	go.draw()
	go.moveTo(Point(1,1))
	go.erase()
