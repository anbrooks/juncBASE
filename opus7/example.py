#!/usr/bin/python2.4
#
#   Copyright (c) 2003 by Bruno R. Preiss, P.Eng.
#
#   $Author: brpreiss $
#   $Date: 2005/06/09 00:00:38 $
#   $RCSfile: example.py,v $
#   $Revision: 1.12 $
#
#   $Id: example.py,v 1.12 2005/06/09 00:00:38 brpreiss Exp $
#

"""
Provides a number of example functions.
"""

__author__  = "Bruno R. Preiss, P.Eng."
__date__    = "$Date: 2005/06/09 00:00:38 $"
__version__ = "$Revision: 1.12 $"
__credits__ = "Copyright (c) 2003 by Bruno R. Preiss, P.Eng."

import sys
import math
import exceptions
from opus7.denseMatrix import DenseMatrix
from opus7.randomNumberGenerator import RandomNumberGenerator

#{
def sum(n):
    result = 0
    i = 1
    while i <= n:
        result += i
        i += 1
    return result
#}>a

#{
#!def Horner(a, n, x):
#[
def Horner1(a, n, x):
#]
    result = a[n]
    i = n - 1
    while i >= 0:
        result = result * x + a[i]
        i -= 1
    return result
#}>b

#{
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)
#}>c

#{
def findMaximum(a, n):
    result = a[0]
    i = 1
    while i < n:
        if a[i] > result:
            result = a[i]
        i += 1
    return result
#}>d

#{
def gamma():
    result = 0.
    i = 1
    while i <= 500000:
        result += 1.0/i - math.log((i + 1.0)/i)
        i += 1
    return result
#}>e

#{
#!def geometricSeriesSum(x, n):
#[
def geometricSeriesSum1(x, n):
#]
    sum = 0
    i = 0
    while i <= n:
        prod = 1
        j = 0
        while j < i:
            prod *= x
            j += 1
        sum += prod
        i += 1
    return sum
#}>f

#{
#!def geometricSeriesSum(x, n):
#[
def geometricSeriesSum2(x, n):
#]
    sum = 0
    i = 0
    while i <= n:
        sum = sum * x + 1
        i += 1
    return sum
#}>g

#{
def power(x, n):
    if n == 0:
        return 1
    elif n % 2 == 0: # n is even
        return power(x * x, n / 2)
    else: # n is odd
        return x * power(x * x, n / 2)
#}>h

#{
#!def geometricSeriesSum(x, n):
#[
def geometricSeriesSum3(x, n):
#]
    return (power(x, n + 1) - 1) / (x - 1)
#}>i

#{
#!def Horner(a, n, x):
#[
def Horner2(a, n, x):
#]
    result = a[n]
    i = n - 1
    while i >= 0:
        result = result * x + a[i]
        i -= 1
    return result
#}>j

#{
def prefixSums(a, n):
    j = n - 1
    while j >= 0:
        sum = 0
        i = 0
        while i <= j:
            sum += a[i]
            i += 1
        a[j] = sum
        j -= 1
#}>k

#{
#!def Fibonacci(n):
#[
def Fibonacci1(n):
#]
    previous = -1
    result = 1
    i = 0
    while i <= n:
        sum = result + previous
        previous = result
        result = sum
        i += 1
    return result
#}>l

#{
#!def Fibonacci(n):
#[
def Fibonacci2(n):
#]
    if n == 0 or n == 1:
        return n
    else:
        #!return Fibonacci(n - 1) + Fibonacci(n - 2)
#[
        return Fibonacci2(n - 1) + Fibonacci2(n - 2)
#]
#}>m

#{
def bucketSort(a, n, buckets, m):
    for j in range(m):
        buckets[j] = 0
    for i in range(n):
        buckets[a[i]] += 1
    i = 0
    for j in range(m):
        for k in range(buckets[j]):
            a[i] = j
            i += 1
#}>n

#{
def binarySearch(array, target, i, n):
    if n == 0:
        raise KeyError
    if n == 1:
        if array[i] == target:
            return i
        raise KeyError
    else:
        j = i + n / 2
        if array[j] <= target:
            return binarySearch(array, target, j, n - n/2)
        else:
            return binarySearch(array, target, i, n/2)
#}>o

#{
#!def Fibonacci(n):
#[
def Fibonacci3(n):
#]
    if n == 0 or n == 1:
        return n
    else:
        #!a = Fibonacci((n + 1) / 2)
        #!b = Fibonacci((n + 1) / 2 - 1)
#[
        a = Fibonacci3((n + 1) / 2)
        b = Fibonacci3((n + 1) / 2 - 1)
#]
        if n % 2 == 0:
            return a * (a + 2 * b)
        else:
            return a * a + b * b
#}>p

def merge(array, pos, m, n):
    temp = [ None ] * (m + n)
    i = pos
    left = pos + m
    j = left
    right = left + n
    k = 0
    while i < left and j < right:
        if array[i] < array[j]:
            temp[k] = array[i]
            k += 1
            i += 1
        else:
            temp[k] = array[j]
            k += 1
            j += 1
    while i < left:
        temp[k] = array[i]
        k += 1
        i += 1
    while j < right:
        temp[k] = array[j]
        k += 1
        j += 1
    for k in xrange(m + n):
        array[pos + k] = temp[k]

#{
def mergeSort(array, i, n):
    if n > 1:
        mergeSort(array, i, n / 2)
        mergeSort(array, i + n / 2, n - n / 2)
        merge(array, i, n / 2, n - n / 2)
#}>q

#{
#!def Fibonacci(n, k):
#[
def Fibonacci4(n, k):
#]
    if n < k - 1:
        return 0
    elif n == k - 1:
        return 1
    else:
        f = [0] * (n + 1)
        for i in xrange(k - 1):
            f[i] = 0
        f[k - 1] = 1
        for i in xrange(k, n + 1):
            sum = 0
            for j in xrange(k + 1):
                sum += f[i - j]
            f[i] = sum
        return f[n]
#}>r

#{
def binom(n, m):
    b = [0] * (n + 1)
    b[0] = 1
    for i in xrange(1, n + 1):
        b[i] = 1
        j = i - 1
        while j > 0:
            b[j] += b[j - 1]
            j -= 1
    return b[m]
#}>s

#{
def typeset(l, D, s):
    n = len(l)
    L = DenseMatrix(n, n)
    for i in xrange(n):
        L[i, i] = l[i]
        for j in xrange(i + 1, n):
            L[i, j] = L[i, j - 1] + l[j]
    P = DenseMatrix(n, n)
    for i in xrange(n):
        for j in xrange(i, n):
            if L[i, j] < D:
                P[i, j] = abs(D - L[i, j] - (j - i) * s)
            else:
                P[i, j] = sys.maxint
    c = DenseMatrix(n, n)
    for j in xrange(n):
        c[j, j] = P[j, j]
        i = j - 1
        while i >= 0:
            min = P[i, j]
            for k in xrange(i, j):
                tmp = P[i, k] + c[k + 1, j]
                if tmp < min:
                    min = tmp
            c[i, j] = min
            i -= 1
#[    
    for i in xrange(n):
        for j in xrange(i, n):
            print c[i, j],
        print
#]
#}>t

#{
def pi(trials):
    hits = 0
    for i in xrange(trials):
        x = RandomNumberGenerator.next
        y = RandomNumberGenerator.next
        if x * x + y * y < 1.0:
            hits += 1
    return 4.0 * hits / trials
#}>u

#{
def one():
    x = 1
    print x
    two(x)
    print x

def two(y):
    print y
    y = 2
    print y
#}>v

#{
class A(Exception):
    pass

def f():
    raise A

def g():
    try:
        f()
    except A:
        # ...
#[
        print "caught A"
#]
#}>w

if __name__ == "__main__":
    print "sum(10) = ", sum(10)
    print "Horner1([2,4,6], 2, 57) = ", Horner1([2,4,6], 2, 57)
    print "Horner2([2,4,6], 2, 57) = ", Horner2([2,4,6], 2, 57)
    print "factorial(10) =", factorial(10)
    print "findMaximum([3,1,4,1,5,9,2], 7) = ", findMaximum([3,1,4,1,5,9,2], 7)
    print "gamma() = ", gamma()
    print "geometricSeriesSum1(10, 6) = ", geometricSeriesSum1(10, 6)
    print "geometricSeriesSum2(10, 6) = ", geometricSeriesSum2(10, 6)
    print "geometricSeriesSum3(10, 6) = ", geometricSeriesSum3(10, 6)
    arg = [2,4,6,8]
    prefixSums(arg, 4)
    print "prefixSums([2,4,6,8], 4) = ", arg
    print "Fibonacci1(5) =", Fibonacci1(10)
    print "Fibonacci2(5) =", Fibonacci2(10)
    print "Fibonacci3(5) =", Fibonacci3(10)
    print "Fibonacci4(5, 2) =", Fibonacci4(10, 2)
    arg = [3,1,4,1,5,9,2]
    bucketSort(arg, len(arg), [0] * 10, 10)
    print "bucketSort([3,1,4,1,5,9,2], 10) = ", arg
    arg = [3,1,4,1,5,9,2]
    mergeSort(arg, 0, 7)
    print "mergeSort([3,1,4,1,5,9,2], 10) = ", arg
    print "binarySearch([1,1,2,3,4,5,9], 5, 0, 7) = ", \
        binarySearch([1,1,2,3,4,5,9], 5, 0, 7)
    print "binom(5, 2) = ", binom(5, 2)
    print "pi(10000) = ", pi(10000)
    one()
