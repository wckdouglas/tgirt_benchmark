#/usr/bin/env python

import pandas as pd
import os

project_path = os.environ['WORK'] + '/cdw2854/bench_marking_new/bench_marking/genome_mapping'
count_path = project_path + '/conventional/counts'

df = pd.read_table(count_path + '/counts.tsv',
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
df.to_csv(tablename, index=False, sep='\t')
print('Written %s' %tablename)

