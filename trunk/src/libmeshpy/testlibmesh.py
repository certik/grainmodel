#! /usr/bin/env python

from distutils.util import get_platform
import os
import sys

import numpy

libdir=os.path.join("build", "lib.%s-%s" % (get_platform(), sys.version[:3]) )
sys.path.append(libdir)
import libmeshpy
sys.path.append(os.path.join("/home/ondra/libmeshpetscpackage/libs/petsc4py",
    libdir))
import petsc4py
petsc4py.init(sys.argv,"linux")
from petsc4py.PETSc import Mat, KSP


import progressbar

print libmeshpy.doubleSum([1,2,3.14,4]), 10.14

myArray = numpy.zeros(5,'d')
libmeshpy.doubleOnes(myArray)
print myArray, numpy.array([1.,1.,1.,1.,1.])

print "starting"
#libmeshpy.mesh("../../tmp/in.xda")

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
    def loadsize(self):
        return self.m.readint()

class system:
    def load(self,fname):
        l=load(fname)
        n=l.loadsize()
        widgets=['Assembling: ', progressbar.Percentage(), ' ', 
                progressbar.Bar(), ' ', progressbar.ETA()]
        self.pbar=progressbar.ProgressBar(widgets=widgets,maxval=n).start()

        self.A=Mat()
        self.A.create()
        self.A.setSizes(n)
        self.A.setFromOptions()
        self.x,self.b=self.A.getVecs()
        for i in range(n):
            Ae=l.loadmatrix()
            l.loadindices()
            Fe=l.loadvector()
            indices=l.loadindices()
            self.A.setValues(indices,indices,Ae)
            self.b.setValues(indices,Fe)
            self.pbar.update(i)

        self.A.assemble()
        #temporary
#        self.b.setRandom()
#        self.A.diagonalSet(self.b)
    def solve(self):
        ksp = KSP()
        ksp.create()
        ksp.setOperators(self.A,self.A,Mat.Structure.SAME_NONZERO_PATTERN)
        ksp.setFromOptions()
        ksp.solve(self.b,self.x)

s=system()
s.load("../../tmp/matrices")
print "solve"
s.solve()
#print s.x.view()

print "ok, we are done."
