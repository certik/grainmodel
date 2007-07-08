#! /usr/bin/env python

from distutils.core import setup, Extension
import numpy
from sys import argv

if argv[1] == "build":
    argv[1] = "build_ext"
if argv[1] == "build_ext":
    argv.insert(2, "--swig-opts=-c++")

extdir="femlib/libmeshpy"
_libmeshpy = Extension("femlib.libmeshpy._libmeshpy",
                    [extdir+"/libmeshpy.i", extdir+"/libmeshpy.cxx"],
                    include_dirs = [numpy.get_include(),"/usr/include/libmesh",
                        "/usr/include/mpi","/usr/include/petsc","femlib/src"],
                    extra_compile_args=["-O2"],
                    libraries = ["mesh","petsc","petscdm","petscksp","petscmat",
                        "petscsnes","petscts","petscvec"],
                    )

setup(name        = "libmeshpy",
      description = "Libmesh bindings for python",
      author      = "Ondrej Certik",
      packages = ["femlib", "geom", "femlib.libmeshpy"],
      #py_modules  = ["femlib.libmeshpy"],
      ext_modules = [_libmeshpy],
      )
