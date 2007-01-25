#! /usr/bin/python
import spmatrix, jdsym, itsolvers, precon
import Numeric

def readmat(f,n):
    A=Numeric.zeros((n,n),"f")
    for j in range(n):
        row=[float(x) for x in f.readline().split()]
        for i in range(n):
            A[j,i]=row[i]
    return A

def gettopo(fname):
    topo=[]
    for l in file(fname):
        indices=[int(x) for x in l.split()]
        topo.append(indices[1:])
    return topo

def update_add_mask_sym2(A,B,ind,mask):
    for i in range(len(ind)):
        for j in range(len(ind)):
            if mask[i]:
                A[ind[i],ind[j]]+=B[i,j]

def assemble(fname,topo,mapping):
    nnodes=len(mapping)
    f=file(fname)
    A=spmatrix.ll_mat(nnodes,nnodes)
    M=spmatrix.ll_mat(nnodes,nnodes)
    for el in topo:
        n=len(el)
        matK=readmat(f,n)
        matV=readmat(f,n)
        matM=readmat(f,n)
        el=[mapping[x] for x in el]
        A.update_add_mask_sym(matK,Numeric.array(el),Numeric.ones(n))
        A.update_add_mask_sym(matV,Numeric.array(el),Numeric.ones(n))
        M.update_add_mask_sym(matM,Numeric.array(el),Numeric.ones(n))
    return A,M

def setzeronodes(A,M,fname,mapping,lam=10000):
    nodes=[int(x)-1 for x in file(fname).readline().split()]
    nodes=[mapping[x] for x in nodes]
    penalty=1e10
    for nod in nodes:
        A[nod,nod]=lam*penalty
        M[nod,nod]=penalty

def remap(vec,map):
#    return [vec[map[x]] for x in range(len(map))]
    return [vec[x] for x in map]

def savevec(vec,fname):
    f=file(fname,"w")
    for n in vec:
        f.write("%f "%n)
    f.write("\n")

print "loading..."
mapping=[int(x) for x in file("tmp/nodemap.libmesh").readline().split()]
#mapping=range(len(mapping))
topo=gettopo("tmp/topo.dat")
A,M=assemble("tmp/matrices.dat",topo,mapping)
setzeronodes(A,M,"tmp/zeronodes.gmsh",mapping)
print "saving..."
A.export_mtx("tmp/matA.dat")
M.export_mtx("tmp/matM.dat")
print "solving..."
tau=0.0
#tau=-10
Atau=A.copy()
Atau.shift(-tau,M)
K=precon.jacobi(Atau)
nE=100
#K=None
A=A.to_sss();M=M.to_sss();
kconv, lmbd, Q, it, it_in = jdsym.jdsym(A, M, K, nE, tau, 1e-5, 150, 
        itsolvers.qmrs, clvl=1, strategy=1)

print "number of converged eigenvalues:",kconv
print lmbd
assert len(Q[:,0])==len(mapping)
for i in range(kconv):
    savevec(remap(Q[:,i],mapping),"out/sol-%d.dat"%i)
    #savevec(Q[:,i],"out/sol-%d.dat"%i)

print "done"
