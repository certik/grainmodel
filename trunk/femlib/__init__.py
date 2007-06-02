from distutils.util import get_platform
import os
import sys

libdir=os.path.join("build", "lib.%s-%s" % (get_platform(), sys.version[:3]) )
sys.path.append(os.path.join("/home/ondra/grainmodel/femlib/src",libdir))

from fem import System, EM
