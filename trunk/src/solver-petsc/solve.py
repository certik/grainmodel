#! /usr/bin/python

import sys, os

if len(sys.argv) == 1:
    #options="-eps_type arpack -eps_smallest_real -eps_nev 2 -eps_tol 1e-3 -eps_view"
    options="-st_type sinvert"
else:
    options=" ".join(sys.argv[1:])

os.system("src/solver-petsc/solver -f1 ../../tmp/matA.matlab -f2 ../../tmp/matM.matlab %s"%(options))
