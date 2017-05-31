#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import cjson
import os
import re

def read_sample(sample_folder):
    samplename = os.path.basename(sample_folder)
    json_file = open(sample_folder + '/lib_format_counts.json')
    salmon = cjson.decode(json_file.read())
    df = pd.DataFrame({'variable':salmon.keys(),'value':salmon.values() }) \
            .drop(6) \
            .assign(sample = samplename)
    return df

project_path = '/stor/work/Lambowitz/cdw2854/bench_marking/'
datapath = project_path + '/alignment_free/salmon'
samples = glob.glob(datapath + '/Sample-*')
samples = filter(lambda x: re.search('-[123]$',x), samples)
df = pd.concat(map(read_sample, samples), axis=0) \
    .query("value!=0") \
    .pipe(pd.pivot_table, index='sample', 
          columns = 'variable', values='value', 
          aggfunc=np.sum)  \
    .reset_index()

tablename = project_path + '/tables/salmon_table.tsv'
df.to_csv(tablename, sep='\t', index=False)
print 'Written %s' %tablename

