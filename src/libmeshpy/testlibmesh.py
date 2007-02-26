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
from petsc4py.PETSc import Mat, KSP, InsertMode


import progressbar

#print libmeshpy.doubleSum([1,2,3.14,4]), 10.14

#myArray = numpy.zeros(5,'d')
#libmeshpy.doubleOnes(myArray)
#print myArray, numpy.array([1.,1.,1.,1.,1.])

print "starting"
libmeshpy.mesh("../../tmp/in.xda")

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
        return self.m.readint(),self.m.readint()

class system:
    def load(self,fname):
        l=load(fname)
        nn,ne=l.loadsize()
        print "nodes:",nn
        print "elements:",ne
        self.nele=ne
        widgets=['Assembling: ', progressbar.Percentage(), ' ', 
                progressbar.Bar(), ' ', progressbar.ETA()]
        self.pbar=progressbar.ProgressBar(widgets=widgets,maxval=ne-1).start()

        IM=InsertMode.ADD_VALUES

        self.A=Mat()
        self.A.create()
        self.A.setSizes(nn)
        self.A.setFromOptions()
        self.x,self.b=self.A.getVecs()
        self.b.zeroEntries()
        for i in range(ne):
            indices=l.loadindices()
            Ae=l.loadmatrix()
            Fe=l.loadvector()
            self.A.setValues(indices,indices,Ae,IM)
            self.b.setValues(indices,Fe,IM)
            self.pbar.update(i)
        self.A.assemble()

    def solve(self):
        ksp = KSP()
        ksp.create()
        ksp.setOperators(self.A,self.A,Mat.Structure.SAME_NONZERO_PATTERN)
        ksp.setFromOptions()
        ksp.solve(self.b,self.x)

s=system()
s.load("../../tmp/matrices")
print "solving"
s.solve()

print "gradient"
x = s.x.getArray()
g = numpy.zeros(s.nele,'d')
libmeshpy.grad("../../tmp/in.xda",x,g)

print "integrating"
print "integ=", libmeshpy.integ("../../tmp/in.xda",g)

print "saving"
import tables
h5=tables.openFile("../../tmp/sol.h5",mode="w",title="Test")
gsol=h5.createGroup(h5.root,"solver","Ax=b")
h5.createArray(gsol,"x",x,"solution vector")
h5.createArray(gsol,"grad",g,"gradient of the solution vector")
h5.close()


print "ok, we are done."
