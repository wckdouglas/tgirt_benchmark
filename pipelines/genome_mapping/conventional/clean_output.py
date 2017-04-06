#/usr/bin/env python


import pandas as pd
import os

count_path = '/scratch/02727/cdw2854/bench_marking/genome_mapping/pipeline7/conventional'

df = pd.read_table(count_path + '/counts.tsv',
        skiprows=1) 
df.columns = map(lambda x: os.path.basename(x.replace('/Hisat/hisat.bam','')), df.columns)
df.drop(['Chr','Start','End','Strand','Length'], axis=1, inplace=True)
df.rename(columns = {'Geneid':'id'},inplace=True)
df.to_csv(count_path + '/feature_counts.tsv', index=False, sep='\t')

