#!/usr/bin/env python

from Bio import Entrez, SeqIO
import sys
from itertools import izip

out_prefix = sys.argv[1]
rRNA_gene = ['gi|23898|emb|X12811.1| Human 5S DNA',
    'gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit']
Entrez.email = 'wckdouglas@gmail.com'

rRNA_fa = out_prefix + '.fa'
rRNA_bed = out_prefix + '.bed'
with open(rRNA_fa,'w') as fa:
    for rRNA in rRNA_gene:
        id = rRNA.split('|')[1]
        handle = Entrez.efetch(db="nucleotide", id=id, 
                           rettype="fasta", retmode="text")
        record = handle.read()
        fa.write('>' + rRNA + '\n'.join(record.split('\n')[1:]) + '\n')

with  open(rRNA_bed,'w') as bed:
    names = '5S rRNA','18S/28S/5.8S rRNA'
    for record,name in izip(SeqIO.parse(rRNA_fa,'fasta'),names):
        bed.write('\t'.join([record.id, '0', str(len(record.seq)), name,'0','+','rRNA'])+'\n')
print 'Written: %s and %s' %(rRNA_fa, rRNA_bed)
    
