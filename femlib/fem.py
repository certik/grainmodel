#! /usr/bin/env python

import sys

import numpy

import libmeshpy
import petsc4py
petsc4py.init(sys.argv,"linux")
from petsc4py.PETSc import Mat, KSP, InsertMode
import progressbar

class MyBar(libmeshpy.Updater):
    """Encapsulation of a nice progress bar"""

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
    """Encapsulates the loadmatrices C++ class for loading the element
    matrices stored in a binary format"""

    def __init__(self,fname):
        self.m=libmeshpy.loadmatrices(fname)

    def loadmatrix(self):
        "Reads 16 floats and returns a 4x4 numpy array"
        x=numpy.zeros((4,4))
        for j in range(4):
            for i in range(4):
                x[i,j]=self.m.readfloat()
        return x

    def loadindices(self):
        "Reads (int,int,int,int) and returns a numpy vector"
        x=numpy.zeros((4),dtype=int)
        for i in range(4):
            x[i]=self.m.readint()
        return x

    def loadvector(self):
        "Reads (float,float,float,float) and returns a numpy vector"
        x=numpy.zeros((4))
        for i in range(4):
            x[i]=self.m.readfloat()
        return x

    def loadsize(self):
        "Reads (int,int)  = (number of nodes, number of elements)"
        return self.m.readint(),self.m.readint()

class System:

    def __init__(self,fmesh,fmatrices,fboundaries):
        self.fmesh=fmesh
        self.fmatrices=fmatrices
        self.fboundaries=fboundaries

    def compute_element_matrices(self):
        """Calculates the element matrices and RHS and includes boundary
        conditions in them. The matrices and RHS together with global indices
        and are stored in self.fmatrices file."""
        libmeshpy.mesh(self.fmesh,self.fmatrices,self.fboundaries,
                MyBar("Element matrices and RHS: "))

    def assemble(self):
        """Ax=b    Assembles A and b using the matrices and RHS and indices
        stored in the file self.fmatrices.
        
        self.A, self.b and self.x then contain initialized petsc matricx and
        vectors.
        """
        l=Matrices(self.fmatrices)
        nn,ne=l.loadsize()
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

        The self.x is also returned from "solve", so that the user doesn't have
        to know about self.x.
        
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
        return self.x

    def gradient(self,x):
        """Computes a gradient of the scalar field defined on every
        node (numpy array x). 
        
        The gradient is computed on every element and the result is returned as
        a list of 3 numpy arrays (x,y,z) components."""
        gx = numpy.zeros(self.nele,'d')
        gy = numpy.zeros(self.nele,'d')
        gz = numpy.zeros(self.nele,'d')
        libmeshpy.grad(self.fmesh,x,gx,gy,gz,MyBar("Gradient: "))
        return gx,gy,gz

    def integ(self, grad, boundarynum):
        """Integrates the normal component of "grad" (list of 3 numpy arrays of
        floats for each element) over the boundary given by "boundarynum".
        
        Returns a float (=the value of the surface integral)."""
        x,y,z=grad
        return libmeshpy.integ(self.fmesh,self.fboundaries,x,y,z,boundarynum,
            MyBar("Integrating over surface %d: "%(boundarynum)))

    def save(self, fsol, x, g):
        "Saves x and g to a file fsol"
        import tables
        h5=tables.openFile(fsol,mode="w",title="Test")
        gsol=h5.createGroup(h5.root,"solver","Ax=b")
        h5.createArray(gsol,"x",x,"solution vector")
        h5.createArray(gsol,"grad",g,"gradient of the solution vector")
        h5.close()
