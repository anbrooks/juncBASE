#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:41 $
#   $RCSfile: timer.py,v $
#   $Revision: 1.6 $
#
#   $Id: timer.py,v 1.6 2005/06/09 00:00:41 brpreiss Exp $
#

"""
Provides a class for measuring execution time.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:41 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import time
from opus7.exception import *

class Timer(object):
    """
    A timer for measuring execution time.
    """

    STOPPED = 1

    RUNNING = 2

    TOLERANCE = 100

    def __init__(self):
        """
        (Timer) -> None

        Constructs a timer.
        """
        self._startTime = time.time()
        self._stopTime = self._startTime
        self._startClock = time.clock()
        self._stopClock = self._startClock
        self._state = self.STOPPED

    def start(self):
        """
        (Timer) -> None

        Starts this timer.
        """
        if self._state != self.STOPPED:
            raise StateError
        self._startTime = time.time()
        self._startClock = time.clock()
        self._state = self.RUNNING

    def stop(self):
        """
        (Timer) -> None

        Stops this timer.
        """
        if self._state != self.RUNNING:
            raise StateError
        self._stopClock = time.clock()
        self._stopTime = time.time()
        self._state = self.STOPPED

    def getElapsedTime(self):
        """
        (Timer) -> double

        Returns the elapsed time.
        """
        if self._state == self.RUNNING:
            self._stopClock = time.clock()
            self._stopTime = time.time()
        elapsedTime = self._stopTime - self._startTime
        # This is a gross kludge.
        # The clock method measures CPU time
        # but the time method measures real-time
        # so these times are not really interchangeable.
        if elapsedTime < self.TOLERANCE:
            elapsedTime = self._stopClock - self._startClock
        return elapsedTime

    @staticmethod
    def main(*argv):
        "Timer test program."
        print Timer.main.__doc__
        t = Timer()
        t.start()
        for i in range(4000000):
            pass
        t.stop()
        print "Elapsed time %g." % (t.getElapsedTime())
        return 0

if __name__ == "__main__":
    sys.exit(Timer.main(*sys.argv))
