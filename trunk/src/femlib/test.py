"""
Tests that python-petsc and libmespy is working correctly.
"""

from distutils.util import get_platform
import sys,os
libdir=os.path.join("build", "lib.%s-%s" % (get_platform(), sys.version[:3]) )
sys.path.append(os.path.join("/home/ondra/grainmodel/femlib/src",libdir))

import numpy

def petsc():
    print "starting petsc"
    from petsc4py.PETSc import Mat
    A=Mat()
    A.createSeqAIJ(1000)
    print "finishing petsc"

def libmesh():
    print "starting libmesh"
    import libmeshpy
    lam=numpy.zeros((10000),"d")
    libmeshpy.mesh("../tmp/in.xda","../tmp/matrices","../tmp/t12.boundaries",
                [],[],lam,None)
    print "finishing libmesh"

petsc()
petsc()
libmesh()
libmesh()
petsc()
libmesh()
petsc()
