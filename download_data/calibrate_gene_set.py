#!/usr/bin/env python


import pandas as pd
import os, sys

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <ref path>' %sys.argv[0])
    
gene_path = sys.argv[1]

def transcript_table():
    tablename = gene_path + '/transcripts.tsv'
    return pd.read_table(tablename) \
            .pipe(lambda d: d[['gene_id','name','type']]) \
            .drop_duplicates()

def genes_bed():
    tablename = gene_path + '/genes.bed'
    return pd.read_table(tablename,
                         names = ['chrom','start','end','gene_name',
                                  'score','strand','gene_type','gene_id'])\
            .pipe(lambda d: d[['gene_name','gene_type','gene_id']])

    

t_tab = transcript_table()
bed = genes_bed()
df = t_tab.merge(bed, how='inner', on ='gene_id')
union_gene = gene_path + '/union_genes.tsv'
df.to_csv(union_gene, index=False, sep='\t')
print 'Written %s' %union_gene
