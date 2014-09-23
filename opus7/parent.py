#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: parent.py,v $
#   $Revision: 1.5 $
#
#   $Id: parent.py,v 1.5 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the Parent class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.5 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.person import Person

#{
class Parent(Person):
    """
    Represents a parent.
    """

    def __init__(self, name, sex, children):
        """
        (Parent, str, int, (Person, Person, ...)) -> None
        Constructs a parent with the given name and sex
        and collection of children.
        """
        super(Parent, self).__init__(name, sex)
        self._children = children

    def getChild(self, i):
        """
        (Parent, int) -> Person
        Returns the specified child of this parent.
        """
        return self._children[i]

    def __str__(self):
        """
        (Parent) -> str
        Returns a string representation of this parent.
        """
        # ...
#[
        pass
#]
#}>a
