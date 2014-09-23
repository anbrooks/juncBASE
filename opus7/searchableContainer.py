#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: searchableContainer.py,v $
#   $Revision: 1.22 $
#
#   $Id: searchableContainer.py,v 1.22 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SearchableContainer class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.22 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.abstractmethod import abstractmethod
from opus7.container import Container
from opus7.visitor import Visitor

#{
class SearchableContainer(Container):
    """
    Base class from which all searchable container classes are derived.
    """

    def __init__(self):
        """
        (SearchableContainer) -> None
        Constructor.
        """
        super(SearchableContainer, self).__init__()

    @abstractmethod
    def __contains__(self, obj): pass

    @abstractmethod
    def insert(self, obj): pass

    @abstractmethod
    def withdraw(self, obj): pass

    @abstractmethod
    def find(self, obj): pass
#}>a

    @staticmethod
    def main(*argv):
        "SearchableContainer test program."
        print SearchableContainer.main.__doc__
        return 0

if __name__ == "__main__":
    sys.exit(SearchableContainer.main(*sys.argv))
