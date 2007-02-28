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

class MyBar(libmeshpy.Updater):
    def __init__(self,text):
        self.text=text
        libmeshpy.Updater.__init__(self)
    def init(self,max):
        widgets=[self.text, progressbar.Percentage(), ' ', 
                progressbar.Bar(), ' ', progressbar.ETA()]
        self.pbar=progressbar.ProgressBar(widgets=widgets,maxval=max).start()
    def update(self,i):
        self.pbar.update(i)


class Matrices:
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
    def __init__(self,fmesh,fmatrices,fboundaries,fsol):
        self.fmesh=fmesh
        self.fmatrices=fmatrices
        self.fboundaries=fboundaries
        self.fsol=fsol
    def compute_element_matrices(self):
        libmeshpy.mesh(self.fmesh,self.fmatrices,self.fboundaries,
                MyBar("Element matrices and RHS: "))
    def assemble(self):
        l=Matrices(self.fmatrices)
        nn,ne=l.loadsize()
#        print "nodes:",nn
#        print "elements:",ne
        self.nele=ne
        pbar=MyBar("Global matrix and RHS: ")
        pbar.init(ne-1)

        IM=InsertMode.ADD_VALUES

        self.A=Mat()
        self.A.create()
        self.A.setSizes(nn)
        self.A.setFromOptions()
#        self.A.option = Mat.Option.ROWS_SORTED
#        self.A.option = Mat.Option.COLUMNS_SORTED
#        self.A.option = Mat.Option.STRUCTURALLY_SYMMETRIC
        self.x,self.b=self.A.getVecs()
        self.b.zeroEntries()
        for i in range(ne):
            indices=l.loadindices()
            Ae=l.loadmatrix()
            Fe=l.loadvector()
            self.A.setValues(indices,indices,Ae,IM)
            self.b.setValues(indices,Fe,IM)
            pbar.update(i)
        self.A.assemble()

    def solve(self,iterguess=23):
        """Solves the system Ax=b.

        A is stored in self.A
        b is stored in self.b
        x is stored in self.x 

        all three of them must be initialized petsc vectors (and a matrix) and
        the solution will be stored in self.x as a numpy array.
        
        The iterguess is the guess for the number of iterations
        the solver is going to make. This is used in the progress bar. If
        it is overestimated, the progressbar will jump for example from 60% to
        100% at the end, if it is underestimated, the progressbar will stop
        at 99% (exactly at (iterguess-1)/iterguess * 100) but the solver will
        be still computing and the progressbar will be updated to 100% only
        when it returns."""
        pbar=MyBar("Solving Ax=b: ")
        pbar.init(iterguess)
        def upd(ksp,iter,rnorm):
            #print "iter:",iter,"norm:",rnorm
            if iter < iterguess: pbar.update(iter)
        ksp = KSP()
        ksp.create()
        ksp.setOperators(self.A,self.A,Mat.Structure.SAME_NONZERO_PATTERN)
        ksp.setFromOptions()
        ksp.setMonitor(upd)
        ksp.solve(self.b,self.x)
        self.x=self.x.getArray()
        pbar.update(iterguess)

    def gradient(self):
        self.g = numpy.zeros(s.nele,'d')
        libmeshpy.grad(self.fmesh,self.x,self.g,
                MyBar("Gradient of solution: "))

    def integ(self):
        i =[ libmeshpy.integ(self.fmesh,self.fboundaries,self.g,b,
            MyBar("Integrating the gradient (%d): "%(b))) for b in (1,2,3) ]
        print "bottom:", i[1]
        print "top   :", i[0]+i[2],"=",i[0],"+",i[2]

    def save(self):
        import tables
        h5=tables.openFile(self.fsol,mode="w",title="Test")
        gsol=h5.createGroup(h5.root,"solver","Ax=b")
        h5.createArray(gsol,"x",self.x,"solution vector")
        h5.createArray(gsol,"grad",self.g,"gradient of the solution vector")
        h5.close()

s=system("../../tmp/in.xda", "../../tmp/matrices", "../../tmp/t12.boundaries",
        "../../tmp/sol.h5")
s.compute_element_matrices()
s.assemble()
s.solve()
s.gradient()
s.integ()
s.save()

print "done."
