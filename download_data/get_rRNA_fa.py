#!/usr/bin/env python

from __future__ import print_function
from Bio import Entrez, SeqIO
import sys

rRNA_gene = ['gi|23898|emb|X12811.1| Human 5S DNA',
    'gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit']
Entrez.email = 'wckdouglas@gmail.com'

for rRNA in rRNA_gene:
    id = rRNA.split('|')[1]
    handle = Entrez.efetch(db="nucleotide", id=id, 
                           rettype="fasta", retmode="text")
    record = handle.read()
    print('>' + rRNA +'\n'+ '\n'.join(record.split('\n')[1:]))

