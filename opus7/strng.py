#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: strng.py,v $
#   $Revision: 1.19 $
#
#   $Id: strng.py,v 1.19 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides the String class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.19 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import warnings
from opus7.object import Object

# Suppress warnings about left shifts and formating.
warnings.filterwarnings("ignore", "x<<y", FutureWarning)
warnings.filterwarnings("ignore", "%u", FutureWarning)

#{
class String(str, Object):
    """
    String class.
    """

#}@head

#{

    # ...
#}@tail

    def __init__(self, obj):
        """
        (String, object) -> None
        Constructs a string with the string representation of the given object.
        The Object metaclass provides a __new__ method
        that initializes the str instance, so none is defined here.
        """
        pass

    def _compareTo(self, obj):
        """
        (String, String) -> int

        Compares this string with the given string.
        """
        assert isinstance(self, obj.__class__)
        return cmp(str(self), str(obj))

#{
    shift = 6

    mask = ~0 << (31 - shift)

    def __hash__(self):
        """
        (String) -> int
        Hashes this string.
        """
        result = 0
        for c in self:
            result = ((result & String.mask) ^
		result << String.shift ^ ord(c)) & sys.maxint
        return result
#}>a

    @staticmethod
    def testHash():
        "String hash test program."
        print String.testHash.__doc__
        print "ett=0%o" % hash(String("ett"))
        print u"tv\u00e5=0%o" % hash(String("tv\u00e5"))
        print u"tva=0%o" % hash(String("tva"))
        print "tre=0%o" % hash(String("tre"))
        print "fyra=0%o" % hash(String("fyra"))
        print "fem=0%o" % hash(String("fem"))
        print "sex=0%o" % hash(String("sex"))
        print "sju=0%o" % hash(String("sju"))
        print u"\u00e5tta=0%o" % hash(String("\u00e5tta"))
        print u"atta=0%o" % hash(String("atta"))
        print "nio=0%o" % hash(String("nio"))
        print "tio=0%o" % hash(String("tio"))
        print "elva=0%o" % hash(String("elva"))
        print "tolv=0%o" % hash(String("tolv"))
        print "abcdefghijklmnopqrstuvwxy=0%o" % (
            hash(String("abcdefghijklmnopqrstuvwxyz")))
	print "ece.uwaterloo.ca=0%o" % hash(String("ece.uwaterloo.ca"))
	print "cs.uwaterloo.ca=0%o" % hash(String("cs.uwaterloo.ca"))
        print "un=0%o" % hash(String("un"))
        print "deux=0%o" % hash(String("deux"))
        print "trois=0%o" % hash(String("trois"))
        print "quatre=0%o" % hash(String("quatre"))
        print "cinq=0%o" % hash(String("cinq"))
        print "six=0%o" % hash(String("six"))
        print "sept=0%o" % hash(String("sept"))
        print "huit=0%o" % hash(String("huit"))
        print "neuf=0%o" % hash(String("neuf"))
        print "dix=0%o" % hash(String("dix"))
        print "onze=0%o" % hash(String("onze"))
        print "douze=0%o" % hash(String("douze"))

    @staticmethod
    def main(*argv):
        "String test program."
        print String.main.__doc__
        String.testHash()
        return 0

if __name__ == "__main__":
    sys.exit(String.main(*sys.argv))
