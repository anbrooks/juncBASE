#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: hashTable.py,v $
#   $Revision: 1.26 $
#
#   $Id: hashTable.py,v 1.26 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the HashTable class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.26 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.abstractmethod import abstractmethod
from opus7.association import Association
from opus7.strng import String
from opus7.searchableContainer import SearchableContainer

#{
class HashTable(SearchableContainer):
    """
    Base class from which all hash tables are derived.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self):
        """
        (HashTable) -> None
        Constructs this hash table.
        """
        super(HashTable, self).__init__()

    @abstractmethod
    def __len__(self): pass

    def getLoadFactor(self):
        """
        (HashTable) -> double
        Returns the load factor of this hash table.
        """
        return self.count / len(self)

    loadFactor = property(
        fget = lambda self: self.getLoadFactor())
#}>a

#{
    def f(self, obj):
        """
        (HashTable, Object) -> int
        Returns the hash of the given object.
        """
        return hash(obj)

    def g(self, x):
        """
        (HashTable, int) -> int
        Hashes an integer using the division method of hashing.
        """
        return abs(x) % len(self)

    def h(self, obj):
        """
        (HashTable, Object) -> int
        Hashes the specified object
        using the composition of the methods f and g.
        """
        return self.g(self.f(obj))
#}>b

    @staticmethod
    def test(hashTable):
        "HashTable test program."
        print HashTable.test.__doc__
        print hashTable
        hashTable.insert(Association(String("foo"), 12))
        hashTable.insert(Association(String("bar"), 34))
        hashTable.insert(Association(String("foo"), 56))
        print hashTable
        obj = hashTable.find(Association(String("foo")))
        print obj
        hashTable.withdraw(obj)
        print hashTable
        for a in hashTable:
            print a
