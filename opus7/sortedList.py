#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: sortedList.py,v $
#   $Revision: 1.12 $
#
#   $Id: sortedList.py,v 1.12 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the SortedList class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.orderedList import OrderedList

#{
class SortedList(OrderedList):
    """
    Base class from which all sorted list classes are derived.
    """

    def __init__(self):
        """
        (SortedList) -> None
        Constructor.
        """
        super(SortedList, self).__init__()
#}>a

    @staticmethod
    def test(list):
        "SortedList test program."
        print SortedList.test.__doc__
        list.insert(4)
        list.insert(3)
        list.insert(2)
        list.insert(1)
        print list
        obj = list.find(2)
        list.withdraw(obj)
        print list
        for i in list:
            print i
