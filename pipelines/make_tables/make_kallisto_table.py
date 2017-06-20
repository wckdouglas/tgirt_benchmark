#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import cjson
import os

def read_sample(sample_folder):
    samplename = os.path.basename(sample_folder)
    df = pd.read_table(sample_folder + '/abundance.tsv')
    count = df.est_counts.astype('int').sum()
    
    json_file = open(sample_folder + '/run_info.json')
    kallisto = cjson.decode(json_file.read())
    input = kallisto['n_processed']
    return samplename, input, count

project_path = '/stor/work/Lambowitz/cdw2854/bench_marking/'
datapath = project_path + '/alignment_free/kallisto'
samples = glob.glob(datapath + '/Sample-*')
df = map(read_sample, samples)
df = pd.DataFrame(df, columns = ['sample','Trimmed pairs','Features fragments'])
tablename = project_path + '/tables/kallisto_table.csv'
df.columns = df.columns.str.replace('_',' ')
df.sort_values('sample').to_csv(tablename,  index=False)
print 'Written %s' %tablename

