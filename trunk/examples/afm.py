#! /usr/bin/env python

import sys; 
sys.path.append("/home/ondra/femgeom")
sys.path.append("..")

from femlib import EM

a = EM("../tmp")
a.load_geometry("../meshes/afm.geo")
a.meshit()
#a.solve()
