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
        if type == 'rRNA':
            if '5SP' in name:
                name = '5S_rRNA'
            if '5_8S_r' in name or '5-8S' in name:
                name = '5.8S_rRNA'
        print '\t'.join([gene_id, name, tid, type])
