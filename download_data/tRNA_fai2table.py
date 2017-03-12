#!/usr/bin/env python

import re

with open('/stor/work/Lambowitz/ref/GRCh38/tRNA/temp/tRNA.fa','r') as tRNA:
    for line in tRNA:
        if line.startswith('>'):
            seq_id = line.lstrip('>').rstrip('\n')
            seq_name = re.split('[0-9]+',seq_id)[0]
            print '\t'.join([seq_name,seq_name,seq_id,'tRNA'])

