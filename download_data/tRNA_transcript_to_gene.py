#!/usr/bin/env python

from __future__ import print_function
import sys
import re

for line in open(sys.argv[1]):
    if line.startswith('>'):
        tid = line.lstrip('>').rstrip()
        gid = re.sub('[0-9]+-[0-9]$','',tid)
        gid = re.sub('[0-9]+-$','',gid)
        print('%s\t%s\t%s\ttRNA' %(gid, gid, tid))
