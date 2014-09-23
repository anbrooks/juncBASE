#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: point.py,v $
#   $Revision: 1.5 $
#
#   $Id: point.py,v 1.5 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Point class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.5 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

#{
class Point(object):
    """
    Represents a point in in image.
    """

    def __init__(self, x, y):
        """
        (Point, int, int) -> None
        Constructs a point with the given coordinates.
        """
        self._x = x
        self._y = y
    # ...
#}>a

    def __str__(self):
	"""
	(Point) -> str
	Returns a textual representation of this point.
	"""
	return "Point(%d,%d)" % (self._x, self._y)
