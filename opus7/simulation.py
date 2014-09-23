#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:40 $
#   $RCSfile: simulation.py,v $
#   $Revision: 1.12 $
#
#   $Id: simulation.py,v 1.12 2005/06/09 00:00:40 brpreiss Exp $
#

"""
Provides the Simulation class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:40 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.leftistHeap import LeftistHeap
from opus7.association import Association
from opus7.exponentialRV import ExponentialRV

#{
class Simulation(object):
    """
    A discrete-event simulation of an M/M/1 queue.
    """

#}@head

#{

    # ...
#}@tail

#{
    ARRIVAL = 0
    DEPARTURE = 1

    def __init__(self):
        super(Simulation, self).__init__()
        self._eventList = LeftistHeap()
        self._serverBusy = False
        self._numberInQueue = 0
        self._serviceTime = ExponentialRV(100.0)
        self._interArrivalTime = ExponentialRV(100.0)

    def run(self, timeLimit):
        """
        (Simulation, double) -> None
        Runs the simulation up to the given time limit.
        """
        self._eventList.enqueue(self.Event(self.ARRIVAL, 0))
        while not self._eventList.isEmpty:
            evt = self._eventList.dequeueMin()
            t = evt.time
            if t > timeLimit:
                self._eventList.purge()
                break
#[
            print evt
#]
            if evt.type == self.ARRIVAL:
                if not self._serverBusy:
                    self._serverBusy = True
                    self._eventList.enqueue(
                        self.Event(self.DEPARTURE,
                        t + self._serviceTime.next))
                else:
                    self._numberInQueue += 1
                self._eventList.enqueue(self.Event(self.ARRIVAL,
                    t + self._interArrivalTime.next))
            elif evt.type == self.DEPARTURE:
                if self._numberInQueue == 0:
                    self._serverBusy = False
                else:
                    self._numberInQueue -= 1
                    self._eventList.enqueue(
                        self.Event(self.DEPARTURE,
                        t + self._serviceTime.next))
#}>a

#{
    class Event(Association):
        """
        Represents an event in the simulation.
        """

        def __init__(self, type, time):
            """
            (Simulation.Event, int, double):
            Constructs an event of the given type at the given time.
            """
            super(Simulation.Event, self).__init__(time, type)
        
        time = property(
            fget = lambda self: self.getKey())

        type = property(
            fget = lambda self: self.getValue())
#}>b

        def __str__(self):
            """
            (Simulation.Event) -> str
            Returns a string representation of this simulation event.
            """
            if self.type == Simulation.ARRIVAL:
                return "Event {" + str(self.time) + ", arrival}"
            elif self.type == Simulation.DEPARTURE:
                return "Event {" + str(self.time) + ", departure}"
            else:
                return None
