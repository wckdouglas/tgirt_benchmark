#!/bin/env python
from __future__ import print_function
import sys
for line in open(sys.argv[1]):
    if line.startswith('>'):
        tid = line.rstrip('\n').lstrip('>')
        gid = tid.split(':')[0]
        print('%s\t%s\t%s\trRNA' %(gid, gid, tid))
