#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: exception.py,v $
#   $Revision: 1.14 $
#
#   $Id: exception.py,v 1.14 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the InternalError, StateError,
ContainerEmpty and ContainerFull exception classes.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.14 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import exceptions

class InternalError(Exception):
    """
    Raised when situation arises that should never occur.
    """
    pass

class StateError(Exception):
    """
    Raised when an operation on an object is not allowed
    due the state of that object.
    """
    pass

class ContainerEmpty(StateError):
    """
    Raised when a container operation fails because the container is empty.
    """
    pass

class ContainerFull(StateError):
    """
    Raised when a container operation fails because the container is full.
    """
    pass

def main(*argv):
    "Exceptions test program."
    print main.__doc__
    return 0

if __name__ == "__main__":
    sys.exit(main(*sys.argv))
