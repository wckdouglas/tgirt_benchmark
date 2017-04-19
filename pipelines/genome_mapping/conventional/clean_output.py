#/usr/bin/env python

import pandas as pd
import os

project_path = os.environ['WORK'] + '/cdw2854/bench_marking/genome_mapping'
count_path = project_path + '/conventional'

df = pd.read_table(count_path + '/counts.tsv',
        skiprows=1) 
df.columns = map(lambda x: os.path.basename(x.replace('/Hisat/hisat.bam','')), df.columns)
df.columns = map(lambda x: x.replace('-','_'), df.columns)
df.drop(['Chr','Start','End','Strand','Length'], axis=1, inplace=True)
df.rename(columns = {'Geneid':'id'},inplace=True)
colnames = sorted(df.columns)
colnames.remove('id')
new_colnames = ['id']
new_colnames.extend(colnames)
df = df[new_colnames]
df.to_csv(count_path + '/feature_counts.tsv', index=False, sep='\t')

