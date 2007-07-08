#! /usr/bin/env python

import sys

s="".join(open(sys.argv[1]).readlines())
print s%(sys.argv[2])
