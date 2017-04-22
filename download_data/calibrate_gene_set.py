#!/usr/bin/env python


import pandas as pd

gene_path = '/stor/work/Lambowitz/ref/human_transcriptome'

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
df = t_tab.merge(bed, how='outer', on ='gene_id')

