#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: complex.py,v $
#   $Revision: 1.14 $
#
#   $Id: complex.py,v 1.14 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the Complex class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.14 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import math

#{
class Complex(object):
    """
    Sample complex number class.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, real, imag):
        """
        (Complex, double, double) -> None
        Constructs a complex number with the given real and imaginary parts.
        """
        self._real = real
        self._imag = imag
#}>a

#{
    def getReal(self):
        """
        (Complex) -> double
        Returns the real part of this complex number.
        """
        return self._real
    def setReal(self, value):
        """
        (Complex, double) -> double
        Sets the real part of this complex number to the given value
        """
        self._real = value

    real = property(
        fget = getReal,
        fset = setReal)

    def getImag(self):
        """
        (Complex) -> double
        Returns the imaginary part of this complex number.
        """
        return self._imag
    def setImag(self, value):
        """
        (Complex, double) -> double
        Sets the imaginary part of this complex number to the given value
        """
        self._imag = value

    imag = property(
        fget = getImag,
        fset = setImag)
#}>b

#{
    def getR(self):
        """
        (Complex) -> double
        Returns the magnitude of this complex number.
        """
        return math.sqrt(self._real * self._real
	    + self._imag * self._imag)

    def setR(self, value):
        """
        (Complex, double) -> double
        Sets the magnitude of this complex number to the given value.
        """
        theta = self.theta
        self._real = value * math.cos(theta)
        self._imag = value * math.sin(theta)

    r = property(
        fget = getR,
        fset = setR)

    def getTheta(self):
        """
        (Complex) -> double
        Returns the phase of this complex number.
        """
        return math.atan2(self._imag, self._real)

    def setTheta(self, value):
        """
        (Complex, double) -> double
        Sets the phase of this complex number to the given value.
        """
        r = self.r
        self._real = r * math.cos(value)
        self._imag = r * math.sin(value)

    theta = property(
        fget = getTheta,
        fset = setTheta)
#}>c

    def __str__(self):
        """
        (Complex) -> str
        Returns a textual representation of this complex number.
        """
        return "%g+%gi" % (self.real, self.imag)

#{
    def __add__(self, c):
        """
        (Complex, Complex) -> Complex
        Returns a complex number equal to the sum of this complex
        number and the given complex number.
        """
        return Complex(self.real + c.real, self.imag + c.imag)

    def __sub__(self, c):
        """
        (Complex, Complex) -> Complex
        Returns a complex number equal to the difference of this complex
        number and the given complex number.
        """
        return Complex(self.real - c.real, self.imag - c.imag)

    def __mul__(self, c):
        """
        (Complex, Complex) -> Complex
        Returns a complex number equal to the product of this complex
        number and the given complex number.
        """
	return Complex(self.real * c.real - self.imag * c.imag,
	    self.real * c.imag + self.imag * c.real)
#}>d

#{
    @staticmethod
    def main(*argv):
        "Complex test program."
#[
        print Complex.main.__doc__
#]
        c = Complex(0, 0)
	print c
	c.real = 1
	c.imag = 2
	print c
	c.theta = math.pi/2 
	c.r = 50
	print c
	c = Complex(1, 2)
	d = Complex(3, 4)
        print c, d, c+d, c-d, c*d
        return 0
#}>e

if __name__ == "__main__":
    sys.exit(Complex.main(*sys.argv))
