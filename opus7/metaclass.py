#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: metaclass.py,v $
#   $Revision: 1.17 $
#
#   $Id: metaclass.py,v 1.17 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides a metaclass for the Object class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.17 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import sets
import string
from opus7.abstractmethod import abstractmethod

#{
class Metaclass(type):
    """
    Metaclass of the Object class.
    Prevents instantiation of classes that contain abstract methods.
    """

    def __init__(self, name, bases, dict):
        """
        (Metaclass, str, tuple, mapping) -> None

        Initializes this metaclass instance.
        """
        type.__init__(self, name, bases, dict)
        self.__new__ = staticmethod(self.new)

        abstractMethodSet = sets.Set()
        reverseMRO = list(self.__mro__)
        reverseMRO.reverse()
        for cls in reverseMRO:
            for (attrName, attr) in cls.__dict__.iteritems():
                if isinstance(attr, abstractmethod):
                    abstractMethodSet.add(attrName)
                else:
                    abstractMethodSet.discard(attrName)
        self.__abstractmethods__ = list(abstractMethodSet)
        self.__abstractmethods__.sort()

    @staticmethod
    def new(*args, **kwargs):
        """
        (Metaclass, ...) -> object

        Creates an instance of the class using the given arguments.
        Raises a TypeError exception if the class is abstract.
        This method is inserted as the method __new__
        in classes instances derived from Metaclass.
        """
        cls = args[0]
        if len(cls.__abstractmethods__) > 0:
            msg = "Can't instantiate abstract class %s. " % (
                cls.__name__)
            msg += "Missing methods %s." % (
                str(cls.__abstractmethods__))
            raise TypeError, msg
        else:
            for base in cls.__mro__:
                if not isinstance(base, Metaclass) and \
                        base is not type:
                    return base.__new__(*args, **kwargs)
            return object.__new__(*args, **kwargs)
#}>a

    @staticmethod
    def main(*argv):
        "Metaclass test program."
        print Metaclass.main.__doc__

        class Good(object):
            __metaclass__ = Metaclass
        c = Good()

        class String(str, Good):
            pass
        c = String('hello')
        print c

        class Bad( Good):
            def foo(self):
                pass
            foo = abstractmethod(foo)
        try:
            c = Bad()
        except TypeError, msg:
            print "Caught TypeError: %s", str(msg)

        return 0

if __name__ == "__main__":
    sys.exit(Metaclass.main(*sys.argv))
