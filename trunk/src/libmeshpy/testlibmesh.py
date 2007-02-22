#! /usr/bin/env python

from   distutils.util import get_platform
import os
import sys

import numpy as N

libDir = "lib.%s-%s" % (get_platform(), sys.version[:3])
sys.path.insert(0,os.path.join("build", libDir))
import libmeshpy

print libmeshpy.doubleSum([1,2,3.14,4]), 10.14

myArray = N.zeros(5,'d')
libmeshpy.doubleOnes(myArray)
print myArray, N.array([1.,1.,1.,1.,1.])

print "starting"
#libmeshpy.mesh("../../tmp/in.xda")

import numpy
class load:
    def __init__(self,fname):
        self.m=libmeshpy.loadmatrices(fname)
    def loadmatrix(self):
        x=numpy.zeros((4,4))
        for j in range(4):
            for i in range(4):
                x[i,j]=self.m.readfloat()
        return x
    def loadindices(self):
        x=numpy.zeros((4),dtype=int)
        for i in range(4):
            x[i]=self.m.readint()
        return x
    def loadvector(self):
        x=numpy.zeros((4))
        for i in range(4):
            x[i]=self.m.readfloat()
        return x

l=load("../../tmp/matrices")
for i in range(10):
    print l.loadmatrix()
    print l.loadindices()
    print l.loadvector()
    print l.loadindices()

print "ok, we are done."
