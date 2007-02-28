from distutils.util import get_platform
import os
import sys

libdir=os.path.join("build", "lib.%s-%s" % (get_platform(), sys.version[:3]) )
sys.path.append(os.path.join("femlib/src",libdir))
sys.path.append(os.path.join("/home/ondra/libmeshpetscpackage/libs/petsc4py",
    libdir))

from fem import System
