#!/usr/bin/env python

'''
Fix duplicate ID on ensembl transcripts
an example: ENST00000290239.7

For the second time this appears, I will add a .1 at the back of the id
'''

from Bio import SeqIO
import sys

ids = set()
fasta_generator = SeqIO.parse(sys.stdin, 'fasta')
for record in fasta_generator:
    if record.id not in ids:
        id = record.id
    else:
        id = record.id + '.1'
    ids.add(id)
    description = ' '.join(record.description.split(' ')[1:]) 
    print '>%s %s\n%s' %(id, description, record.seq)
