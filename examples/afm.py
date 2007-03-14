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
            #a=0.0001,Q=0.8)
            #a=0.0001,Q=0.8)
            #a=0.00007,Q=0.8,quadratic=True)
            a=0.1,quadratic=True)
    m=geom.read_tetgen("../tmp/t.1")
    m.printinfo()
    m.writemsh("../tmp/t12.msh")
    m.writexda("../tmp/in.xda")
    m.writeBC("../tmp/t12.boundaries")

    global regions,nele
    regions=m.regions
    nele=len(m.elements)

def fem():
    r={100:0.001, 101: 0.1}
    bc={1:0.0, 2: 1.0}

    lam=numpy.zeros((nele),"d")
    for reg in regions:
        for i in regions[reg]:
            lam[i-1]=r[reg]

    s=System("../tmp/in.xda", "../tmp/matrices", "../tmp/t12.boundaries")
    s.compute_element_matrices(bc,lam)
    s.assemble()
    sol=s.solve()
    grad=s.gradient(sol)
    substrate=s.integ(grad,1)
    tip=s.integ(grad,2)
    sides=s.integ(grad,3)
    top=s.integ(grad,4)

    print "results:"
    print "tip      :",tip
    print "substrate:",substrate
    print "top      :",top
    print "sides    :",sides
    print "----- total -----"
    print "all      :",tip+substrate+top+sides

    fname="../tmp/sol.h5"
    m=meshutils.mesh()
    m.readmsh("../tmp/t12.msh")
    m.scalars=sol
    m.writescalarspos(fname[:-4]+".pos","libmesh")
    g=numpy.sqrt(grad[0]**2+grad[1]**2+grad[2]**2)
    m.convert_el_to_nodes(g)
    m.writescalarspos(fname[:-4]+"g.pos","libmesh")

    g = grad[0]**2 + grad[1]**2 + grad[2]**2

    g = convert_grad2constrain(g)

    savevol("../tmp/t.1.vol",g)


    geom.runtetgen("/home/ondra/femgeom/tetgen/tetgen","../tmp/t.1",
            a=0.01,quadratic=True,refine=True)
    m=geom.read_tetgen("../tmp/t.2")
    m.printinfo()
    m.writemsh("../tmp/t12-2.msh")
    #m.writexda("../tmp/in.xda")
    #m.writeBC("../tmp/t12.boundaries")

def savevol(filename,vols):
    f = open(filename,"w")
    f.write("%d\n"%len(vols))
    for i,x in enumerate(vols):
        f.write("%d %f\n"%(i,x))

def convert_grad2constrain(g):
    #max volume constrain
    a = 0.001
    a = 0.1
    g = g/max(g)
    g2 = []
    for x in g:
        if x/a < 1e-9 or a/x > 1.0:
            g2.append(-1)
        else:
            g2.append(a/x)
    return g2

mesh()
fem()
