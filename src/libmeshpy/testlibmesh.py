#! /usr/bin/env python

from   distutils.util import get_platform
import os
import sys

import numpy as N

libDir = "lib.%s-%s" % (get_platform(), sys.version[:3])
sys.path.insert(0,os.path.join("build", libDir))
import libmeshpy

print libmeshpy.doubleSum([1,2,3.14,4]), 10.14

myArray = N.zeros(5,'d')
libmeshpy.doubleOnes(myArray)
print myArray, N.array([1.,1.,1.,1.,1.])

libmeshpy.mesh()
