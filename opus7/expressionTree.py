#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:39 $
#   $RCSfile: expressionTree.py,v $
#   $Revision: 1.6 $
#
#   $Id: expressionTree.py,v 1.6 2005/06/09 00:00:39 brpreiss Exp $
#

"""
Provides the ExpressionTree class.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:39 $"
__version__ = "$Revision: 1.6 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

from opus7.binaryTree import BinaryTree
from opus7.stackAsLinkedList import StackAsLinkedList
from opus7.prePostVisitor import PrePostVisitor

#{
class ExpressionTree(BinaryTree):
    """
    Represents expressions comprised of binary operators.
    """

#}@head

#{

    # ...
#}@tail

#{
    def __init__(self, word):
        """
        Constructs an expression tree with the given word.
        (ExpressionTree, str) -> None
        """
        super(ExpressionTree, self).__init__(word)
    
    @staticmethod
    def parsePostfix(input):
        """
        (File) -> ExpresionTree
        Parses the given file into an expression tree.
        """
        stack = StackAsLinkedList()

        for line in input.readlines():
            for word in line.split():
                if word == "+" or word == "-" \
                        or word == "*" or word == "/":
                    result = ExpressionTree(word)
                    result.attachRight(stack.pop())
                    result.attachLeft(stack.pop())
                    stack.push(result)
                else:
                    stack.push(ExpressionTree(word))
        return stack.pop()
#}>a

#{
    class InfixVisitor(PrePostVisitor):
        """
        Visits the nodes of an expression tree and constructs
        a string that contains the infix representation of the expression.
        """

        def __init__(self):
            """
            (ExpressionTree.InfixVisitor) -> None
            Constructor.
            """
            self._s = ""

        def preVisit(self, obj):
            """
            (ExpressionTree.InfixVisitor, Object) -> None
            Adds a left parenthesis to the string.
            """
            self._s = self._s + "("

        def inVisit(self, obj):
            """
            (ExpressionTree.InfixVisitor, Object) -> None
            Adds the object to the string.
            """
            self._s = self._s + str(obj)

        def postVisit(self, obj):
            """
            (ExpressionTree.InfixVisitor, Object) -> None
            Adds a right parenthesis to the string.
            """
            self._s = self._s + ")"

        def __str__(self):
            """
            (ExpressionTree.InfixVisitor) -> str
            Returns the string.
            """
            return self._s

    def __str__(self):
        """
        (ExpressionTree) -> str
        Returns a string containing the infix representation
        of the expression represented by this expression tree.
        """
        visitor = self.InfixVisitor()
        self.depthFirstTraversal(visitor)
        return str(visitor)
#}>b
