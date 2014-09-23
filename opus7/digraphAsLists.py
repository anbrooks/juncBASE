#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: digraphAsLists.py,v $
#   $Revision: 1.10 $
#
#   $Id: digraphAsLists.py,v 1.10 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides the DigraphAsLists class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.10 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
from opus7.graph import Graph
from opus7.digraph import Digraph
from opus7.graphAsLists import GraphAsLists

class DigraphAsLists(Digraph, GraphAsLists):
    """
    Digraph implemented using adjacency lists.
    """

    def __init__(self, size):
        """
        Constructs a digraph with the given maximum number of vertices.
        """
        super(DigraphAsLists, self).__init__(size)
        self._isDirected = True
    
    @staticmethod
    def main(*argv):
        "DigraphAsLists test program."
        print DigraphAsLists.main.__doc__
        dg = DigraphAsLists(32)
        Graph.test(dg)
        dg.purge()
        Graph.testWeighted(dg)
        return 0

if __name__ == "__main__":
    sys.exit(DigraphAsLists.main(*sys.argv))
