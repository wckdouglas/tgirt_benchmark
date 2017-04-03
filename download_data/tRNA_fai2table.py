#!/usr/bin/env python

import re
import sys

with open(sys.argv[1],'r') as tRNA:
    for line in tRNA:
        if line.startswith('>'):
            seq_id = line.lstrip('>').rstrip('\n')
            seq_name = re.sub('[0-9]+$','',seq_id)
            print '\t'.join([seq_name,seq_name,seq_id,'tRNA'])
