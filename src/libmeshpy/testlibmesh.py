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

print libmeshpy.doubleSum([1,2,3.14,4]), 10.14

myArray = numpy.zeros(5,'d')
libmeshpy.doubleOnes(myArray)
print myArray, numpy.array([1.,1.,1.,1.,1.])

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
        widgets=['Assembling: ', progressbar.Percentage(), ' ', 
                progressbar.Bar(), ' ', progressbar.ETA()]
        self.pbar=progressbar.ProgressBar(widgets=widgets,maxval=ne-1).start()

        #IM=InsertMode.INSERT_VALUES
        IM=InsertMode.ADD_VALUES

        self.A=Mat()
        self.A.create()
        self.A.setSizes(nn)
        self.A.setFromOptions()
        self.x,self.b=self.A.getVecs()
        self.b.zeroEntries()
        for i in range(ne):
            Ae=l.loadmatrix()
            l.loadindices()
            Fe=l.loadvector()
            indices=l.loadindices()
            self.A.setValues(indices,indices,Ae,IM)
#            print indices
            self.b.setValues(indices,Fe,IM)
            self.pbar.update(i)

#        self.A.setValues([5264],[5264],1.0,IM)
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

print "saving"
#f=open("../../tmp/sol.dat","w")
#for a in s.x.getValues(range(len(s.x))):
#    f.write("%f "%a)

import tables
h5=tables.openFile("../../tmp/sol.h5",mode="w",title="Test")
gsol=h5.createGroup(h5.root,"solver","Ax=b")
h5.createArray(gsol,"x",s.x.getArray(),"solution vector")
h5.close()

print "ok, we are done."
