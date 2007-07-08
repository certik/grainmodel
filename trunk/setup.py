#! /usr/bin/env python

from distutils.core import setup, Extension, Command
import numpy
from sys import argv

if argv[1] == "build":
    argv[1] = "build_ext"
if argv[1] == "build_ext":
    argv.insert(2, "--swig-opts=-c++")

class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as 
    is in the svn.
    """
    
    description = "Clean everything"
    user_options = [("all","a","the same")]  

    def initialize_options(self):  
        self.all = None
    
    def finalize_options(self):   
        pass

    def run(self):
        import os
        os.system("py.cleanup")
        os.system("rm -f python-build-stamp-2.4")
        os.system("rm -f MANIFEST")
        os.system("rm -rf build")
        os.system("rm -rf dist")
        os.system("rm -f femlib/libmeshpy/libmeshpy.py")
        os.system("rm -f femlib/libmeshpy/libmeshpy_wrap.cpp")
        os.system("rm -f femlib/libmeshpy/libmeshpy_wrap.h")

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
      cmdclass    = {
                     'clean' : clean, 
                     },
      )
