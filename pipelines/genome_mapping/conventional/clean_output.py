#/usr/bin/env python

import pandas as pd
import numpy as np
import os
import re

def read_genes_info():
    genes_bed = '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/genes.bed'
    return pd.read_table(genes_bed, 
                         usecols = [3,6,7],
                         names = ['name', 'type', 'id'])

def assign_rRNA(row):
    id = row['id']
    if row['type'] == 'rRNA':
        if '5S' in row['name']:
            id = '5S_rRNA' 
        elif '5-8S' in row['name'] or '5.8S' in row['name'] or '5_8S' in row['name']:
            id = '5.8S_rRNA'
        elif '18S' in row['name']:
            id = '18S_rRNA'
        elif '28S' in row['name']:
            id = '28S_rRNA'
    return id



project_path = os.environ['WORK'] + '/cdw2854/bench_marking_new/bench_marking/genome_mapping'
count_path = project_path + '/conventional/counts'

df = pd.read_table(count_path + '/counts',
        skiprows=1) 
df.columns = list(map(lambda x: os.path.basename(x.replace('.bam','')), df.columns))
df.columns = list(map(lambda x: x.replace('-','_'), df.columns))
df.drop(['Chr','Start','End','Strand','Length'], axis=1, inplace=True)
df.rename(columns = {'Geneid':'id'},inplace=True)
colnames = sorted(df.columns)
colnames.remove('id')
new_colnames = ['id']
new_colnames.extend(colnames)
df = df[new_colnames]

tablename =count_path + '/feature_counts.tsv' 
sum_df = df \
    .merge(read_genes_info(), on ='id', how = 'left') \
    .assign(id = lambda d: [assign_rRNA(row) for i, row in d.iterrows()])\
    .assign(id = lambda d: d.id.str.replace('_gene$','')) \
    .assign(id = lambda d: np.where(d.name.str.contains('MT-T'),
                                    d.name.str.replace('[0-9]+$',''),
                                    d.id))\
    .drop(['name','type'], axis=1)\
    .groupby(['id'], as_index=False)\
    .sum() \
    .to_csv(tablename, index=False, sep='\t')
print('Written %s' %tablename)

