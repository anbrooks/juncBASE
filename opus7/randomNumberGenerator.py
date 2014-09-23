#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: randomNumberGenerator.py,v $
#   $Revision: 1.16 $
#
#   $Id: randomNumberGenerator.py,v 1.16 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the RandomNumberGenerator class
and the RandomNumberGenerator singleton.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.16 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys

#{
class RandomNumberGenerator(object):
    """
    A multiplicative linear congruential pseudo-random number generator.

    Adapted from the minimal standard pseudo-random number generator
    described in Stephen K. Park and Keith W. Miller,
    "Random Number Generators: Good Ones Are Hard To Find,"
    Communications of the ACM, Vol. 31, No. 10, Oct. 1988, pp. 1192-1201.
    """

#}@head

#{
#}@tail

#{
    a = 16807
    m = 2147483647
    q = 127773
    r = 2836

    def __init__(self, seed=1):
        """
        (RandomNumberGenerator [, int]) -> None
        Constructs a random number generator with the given seed.
        """
        super(RandomNumberGenerator, self).__init__()
	assert seed >= 1 and seed < self.m
        self._seed = seed

    def getSeed(self):
        """
        (RandomNumberGenerator) -> int
        Returns the seed of this random number generator.
        """
        return self._seed

    def setSeed(self, seed):
        """
        (RandomNumberGenerator, int ) -> None
        Sets the seed of this random number generator to the given value
        """
	assert seed >= 1 and seed < self.m
        self._seed = seed

    seed = property(
        fget = lambda self: self.getSeed(),
        fset = lambda self, value: self.setSeed(value))

    def getNext(self):
        """
        (RandomNumberGenerator) -> double
        Returns the next random number.
        """
        self._seed = self.a * (self._seed % self.q) \
                - self.r * (self._seed / self.q)
        if self._seed < 0:
            self._seed += self.m
        return (1.0 * self._seed)  / self.m

    next = property(
        fget = lambda self: self.getNext())
#[
    @staticmethod
    def main(*argv):
        "RandomNumberGenerator test program."
        print RandomNumberGenerator.main.__doc__
        RandomNumberGenerator.seed = 1
        for i in xrange(10):
            print RandomNumberGenerator.next
        return 0

_RandomNumberGenerator = RandomNumberGenerator
#]

RandomNumberGenerator = RandomNumberGenerator() # singleton
#}>a

if __name__ == "__main__":
    sys.exit(RandomNumberGenerator.main(*sys.argv))
