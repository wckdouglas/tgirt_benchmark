#!/usr/bin/env python

'''
convert ensemble fasta headers to transcript table
'''

import sys
import re
import fileinput

print 'gene_id\tname\tt_id\ttype'
for line in fileinput.input():
    if line.startswith('>'):
        fields = line.lstrip('>').rstrip('\n').split(' ')
        tid = fields[0]
        for field in fields:
            if field.startswith('gene:'):
                gene_id = field.split(':')[1]
                gene_id = re.sub('\.[0-9]+$','',gene_id)
            elif field.startswith('gene_biotype:'):
                type=field.split(':')[1]
            elif field.startswith('gene_symbol:'):
                name = field.split(':')[1]
        print '\t'.join([gene_id, name, tid, type])
                
    
