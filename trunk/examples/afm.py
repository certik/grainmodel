#! /usr/bin/env python

import sys; 
sys.path.append("/home/ondra/femgeom")
sys.path.append("..")

import numpy
import pexpect

import geom
from geom import meshutils
from femlib import System

def mesh():
    meshgeometry="../meshes/afm.geo"

    pexpect.run("gmsh -0 %s -o ../tmp/x.geo"%(meshgeometry))
    g=geom.read_gmsh("../tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"../tmp/t.poly")
    geom.runtetgen("/home/ondra/femgeom/tetgen/tetgen","../tmp/t.poly",
            a=0.001,Q=0.8)
    m=geom.read_tetgen("../tmp/t.1")
    m.printinfo()
    m.writemsh("../tmp/t12.msh")
    m.writexda("../tmp/in.xda")
    m.writeregions("../tmp/t12.regions")
    m.writeBC("../tmp/t12.boundaries")

def fem():
    s=System("../tmp/in.xda", "../tmp/matrices", "../tmp/t12.boundaries")
    s.compute_element_matrices({1: 0.0, 2: 1.0})
    s.assemble()
    sol=s.solve()
    grad=s.gradient(sol)
    substrate=s.integ(grad,1)
    tip=s.integ(grad,2)

    print "results:"
    print "tip      :",tip
    print "substrate:",substrate
    print "----- total -----"
    print "all      :",tip+substrate

    fname="../tmp/sol.h5"
    m=meshutils.mesh()
    m.readmsh("../tmp/t12.msh")
    m.scalars=sol
    m.writescalarspos(fname[:-4]+".pos","libmesh")
    g=numpy.sqrt(grad[0]**2+grad[1]**2+grad[2]**2)
    m.convert_el_to_nodes(g)
    m.writescalarspos(fname[:-4]+"g.pos","libmesh")


mesh()
fem()
