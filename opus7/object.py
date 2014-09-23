#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: object.py,v $
#   $Revision: 1.28 $
#
#   $Id: object.py,v 1.28 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Object class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.28 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.metaclass import Metaclass

#{
class Object(object):
    """
    Base class from which all objects are derived.
    """

#}@head

#{
    # ...
#}@tail

#{
    def __init__(self):
        """
        (Object) -> None
        Constructor.
        """
        super(Object, self).__init__()

    def __cmp__(self, obj):
        """
        (Object, Object) -> int
        Compares this object with the given object.
        """
        if isinstance(self, obj.__class__):
            return self._compareTo(obj)
        elif isinstance(obj, self.__class__):
            return -obj._compareTo(self)
        else:
            return cmp(self.__class__.__name__,
                obj.__class__.__name__)

    @abstractmethod
    def _compareTo(self, obj): pass
#}>a

#{
    __metaclass__ = Metaclass
#}>b

    @staticmethod
    def main(*argv):
        "Object test program."
        print Object.main.__doc__
        try:
            object = Object()
        except TypeError, msg:
            print "Caught TypeError: %s" % str(msg)
        return 0

if __name__ == "__main__":
    sys.exit(Object.main(*sys.argv))
