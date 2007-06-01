#! /usr/bin/env python

from distutils.core import *
import numpy

_libmeshpy = Extension("_libmeshpy",
                    ["libmeshpy_wrap.cxx", "libmeshpy.cxx"],
                    include_dirs = [numpy.get_include(),"/usr/include/libmesh",
                        "/usr/include/mpi","/usr/include/petsc"],
                    extra_compile_args=["-O2"],
                    libraries = ["mesh","petsc","petscdm","petscksp","petscmat",
                        "petscsnes","petscts","petscvec"]
                    )

setup(name        = "libmeshpy",
      description = "Libmesh bindings for python",
      author      = "Ondrej Certik",
      py_modules  = ["libmeshpy"],
      ext_modules = [_libmeshpy]
      )
