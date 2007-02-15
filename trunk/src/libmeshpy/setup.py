#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig
import glob
import os

# Third-party modules - we depend on numpy for everything
import numpy

# _Series extension module
lib="/home/ondra/libmeshpetscpackage/libs/libmesh2"
libmesh=os.path.dirname(glob.glob(lib+"/lib/*/libmesh.so")[0])
libmesh_contrib=os.path.dirname(glob.glob(lib+
    "/contrib/lib/*/liblaspack.so")[0])
libmesh_inc=glob.glob(lib+"/include/*")
libpaths=[libmesh,libmesh_contrib]
_Series = Extension("_libmeshpy",
                    ["libmeshpy_wrap.cxx",
                     "libmeshpy.cxx"],
                    include_dirs = [numpy.get_include()]+libmesh_inc,
                    library_dirs=libpaths,
                    runtime_library_dirs=libpaths,
                    libraries = ["mesh","laspack"]
                    )

# Series setup
setup(name        = "libmeshpy",
      description = "Functions that work on series",
      author      = "Bill Spotz",
      py_modules  = ["libmeshpy"],
      ext_modules = [_Series]
      )
