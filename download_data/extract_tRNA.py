#!/usr/bin/env python

import pandas as pd
import os

ref_path = os.environ['REF'] + '/human_transcriptome'
xlsx_file = ref_path + '/tRNA.xlsx'

tRNA_bed = pd.read_excel(xlsx_file,sheetname=3, header=None) 
bed_file = ref_path + '/tRNA_xlsx.bed'
tRNA_bed.to_csv(bed_file, sep='\t', index=False, header=False)
print 'Written %s' %bed_file



df = pd.read_excel(xlsx_file,sheetname=7,header=None, skiprows=1)   
for i, record in df.iterrows():
    id = record[0]
    if not pd.isnull(id):
        seq = [base for base in record[3:] if not pd.isnull(base)]
        seq = ''.join(seq)
        print '>%s\n%sCCAA' %(id,seq.upper())

mt_tRNA = pd.read_excel(xlsx_file, sheetname=5, header=None)
for i, record in df.iterrows():
    id = record[0]
    seq = record[3]
    print '>%s\n%sCCAA' %(id,seq.upper())

