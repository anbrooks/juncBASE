#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:37 $
#   $RCSfile: abstractmethod.py,v $
#   $Revision: 1.21 $
#
#   $Id: abstractmethod.py,v 1.21 2005/06/09 00:00:37 brpreiss Exp $
#

"""
Provides the abstractmethod descriptor.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:37 $"
__version__ = "$Revision: 1.21 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import inspect
import types

#{
class abstractmethod(object):
    """
    Descriptor for an abstract method.
    """

    def __init__(self, func):
        """
        (abstractmethod, function) -> None

        Constructor.
        """
	assert inspect.isfunction(func)
        self._func = func

    def __get__(self, obj, type):
        """
        (abstractmethod, object, type) -> abstractmethod.method

        Returns an abstractmethod.method that represents
        the binding of this abstract method to the given object.
        """
        return self.method(obj, self._func, type)

    class method(object):
        """
        Abstract method.
        """

        def __init__(self, obj, func, cls):
            """
            (abstractmethod.method, object, function) -> None

            Constructs a binding between the given object instance
            and the given method.
            """
            self._self = obj
            self._func = func
            self._class = cls
            self.__name__ = func.__name__

        def __call__(self, *args, **kwargs):
            """
            (abstractmethod.method, ...) -> None

            Invoked when an abstract method is called.
            """
            msg = "Abstract method %s of class %s called." % (
                    self._func.__name__, self._class.__name__)
            raise TypeError, msg
#}>a

        def __get__(self, *args):
            return self

    @staticmethod
    def main(*argv):
        "abstractmethod test program."
        print abstractmethod.main.__doc__
        
        class Test(object):
            def m(self):
                pass
            m = abstractmethod(m)
            def m2(self, arg):
                pass
            m2 = abstractmethod(m2)

        test = Test()
        try:
            test.m()
        except TypeError, msg:
            print "Caught TypeError: %s" % str(msg)
        try:
            test.m2(57)
        except TypeError, msg:
            print "Caught TypeError: %s" % str(msg)
        print inspect.ismethoddescriptor(test.m)
        print inspect.isroutine(test.m)
        b = test.m.__get__(None, None)
        print inspect.ismethoddescriptor(b)
        print inspect.isroutine(b)
        return 0

if __name__ == "__main__":
    sys.exit(abstractmethod.main(*sys.argv))
