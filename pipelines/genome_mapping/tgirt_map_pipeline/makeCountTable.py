#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import glob
from functools import partial
from multiprocessing import Pool

def readDF(count_file_name):
    df = pd.read_table(count_file_name, header=None)  \
        .pipe(lambda d: d[[3,6,7,8]])
    df.columns = ['name','type','id','count']
    return df

def read_tRNA(count_file_name):
    df = pd.read_table(count_file_name, names= ['id','count']) 
    return df


def readSample(count_file_path, tRNA_count_path, sample_id):    
    print('Running %s' %sample_id)
    df = readDF(count_file_path + '/' + sample_id + '.counts') \
        .pipe(lambda d: d[['id','count']])
    tRNA_df = read_tRNA(tRNA_count_path + '/' + sample_id + '.tRNA')
    df = pd.concat([df, tRNA_df],axis=0) \
        .assign(sample_name = sample_id.replace('-','_')) 
    return df

def main():
    work = os.environ['WORK']
    count_path = work + '/cdw2854/bench_marking/genome_mapping/Counts'
    count_path = work + '/cdw2854/bench_marking_new/bench_marking/genome_mapping/tgirt_map/Counts'
    count_file_path = count_path + '/RAW'
    tRNA_count_path = count_path + '/tRNA_RAW'
    count_files = glob.glob(count_file_path + '/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))
    dfFunc = partial(readSample, count_file_path, tRNA_count_path)
    dfs = Pool(12).map(dfFunc, sample_ids)
    df = pd.concat(dfs, axis=0) \
        .query('sample_name != "try"') \
        .assign(count = lambda d: d['count'].astype(int)) \
        .pipe(pd.pivot_table,index = ['id'],  
            values = 'count' , columns = ['sample_name']) \
        .reset_index()\
        .fillna(0)
    df.iloc[:,1:] = df.iloc[:,1:].astype(int)
    tablename = count_file_path + '/combined_gene_count.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    print('Written %s' %tablename)


if __name__ == '__main__':
    main()
