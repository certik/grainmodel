#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig
import glob

# Third-party modules - we depend on numpy for everything
import numpy

# _Series extension module
lib="/home/ondra/libmeshpetscpackage/libs"
libmeshso=glob.glob(lib+"/libmesh/lib/*/libmesh.so")[0]
petsc="/home/ondra/libmeshpetscpackage/libs/petsc/lib/linux/"
libmesh=libmeshso[:libmeshso.rfind("/")]
libmesh_contrib="/home/ondra/libmeshpetscpackage/libs/libmesh/contrib/lib/i686-pc-linux-gnu_opt"
slepc="/home/ondra/libmeshpetscpackage/libs/slepc/lib/linux/"
_Series = Extension("_Series",
                    ["Series_wrap.cxx",
                     "series.cxx"],
                    include_dirs = [numpy.get_include(),
                        lib+"/libmesh/include/base",
                        lib+"/libmesh/include/utils",
                        lib+"/libmesh/include/enums",
                        lib+"/libmesh/include/mesh",
                        lib+"/petsc/include/mpiuni",
                        lib+"/petsc/bmake/linux/"],
                    library_dirs=[libmesh,petsc,libmesh_contrib,slepc],
                    runtime_library_dirs=[libmesh,petsc,libmesh_contrib,slepc],
                    libraries = ["mesh","petsc","petscts","petscsnes",
                        "petscvec","petscmat","petscdm","petscksp",
                        "laspack","parmetis","metis","sfcurves","gzstream",
                        "tetgen","triangle",
                        "slepc"]
                    )

# Series setup
setup(name        = "Series",
      description = "Functions that work on series",
      author      = "Bill Spotz",
      py_modules  = ["Series"],
      ext_modules = [_Series]
      )
