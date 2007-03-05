#! /usr/bin/env python

import sys; 
sys.path.append("/home/ondra/femgeom")
sys.path.append("..")

import numpy
import pexpect

import geom
from geom import meshutils
from femlib import System

def mesh(tipr=0.5, mesha=0.01, meshQ=0.8, verbose=False):
    s="".join(open("afm-templ.geo").readlines())
    meshgeometry="../tmp/t.geo"
    open(meshgeometry,"w").write(s%(tipr))

    pexpect.run("gmsh -0 %s -o ../tmp/x.geo"%(meshgeometry))
    g=geom.read_gmsh("../tmp/x.geo")
    if verbose: g.printinfo()
    geom.write_tetgen(g,"../tmp/t.poly")
    geom.runtetgen("/home/ondra/femgeom/tetgen/tetgen","../tmp/t.poly",
            a=mesha,Q=meshQ,quadratic=True,verbose=verbose)
    m=geom.read_tetgen("../tmp/t.1",verbose)
    if verbose: m.printinfo()
    m.writemsh("../tmp/t12.msh",verbose)
    m.writexda("../tmp/in.xda",verbose)
    m.writeBC("../tmp/t12.boundaries",verbose)

    regions=m.regions
    nele=len(m.elements)
    return regions, nele

def fem(regions,nele,verbose=False):
    r={100:0.001, 101: 0.1}
    bc={1:0.0, 2: 1.0}

    lam=numpy.zeros((nele),"d")
    for reg in regions:
        for i in regions[reg]:
            lam[i-1]=r[reg]

    s=System("../tmp/in.xda", "../tmp/matrices", "../tmp/t12.boundaries",
            verbose)
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

    #fname="../tmp/sol.h5"
    #m=meshutils.mesh()
    #m.readmsh("../tmp/t12.msh")
    #m.scalars=sol
    #m.writescalarspos(fname[:-4]+".pos","libmesh")
    #g=numpy.sqrt(grad[0]**2+grad[1]**2+grad[2]**2)
    #m.convert_el_to_nodes(g)
    #m.writescalarspos(fname[:-4]+"g.pos","libmesh")
    return substrate

def calc(tipr=0.5,mesha=0.1,meshQ=None,verbose=False):
    r,e=mesh(tipr,mesha,meshQ,verbose)
    return fem(r,e,verbose)

vals=numpy.arange(0.4,0.8,0.05)
#vals=[0.5]
res=[]
for tipr in vals:
    print tipr,"-"*40
    a= calc(tipr,mesha=0.01,meshQ=0.8)
    res+=[a]

print res
