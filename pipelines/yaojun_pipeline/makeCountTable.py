#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import glob
from functools import partial
from multiprocessing import Pool

def readrRNACount(count_file_name):
    rDNA_genes = pd.DataFrame({'id':['gi|23898|emb|X12811.1|','gi|555853|gb|U13369.1|HSU13369'],
                                'type':['rRNA','rRNA'],
                                'name': ['5S_rRNA','rRNA_cluster']})
    return pd.read_table(count_file_name,names=['id','count']) \
            .merge(rDNA_genes, on = 'id', how='inner') 

def readGenes(count_file_name):
    return pd.read_table(count_file_name,names = ['name','type','id','count'])

def readDF(count_file_name):
    if 'rRNA' in count_file_name:
        return readrRNACount(count_file_name)
    elif 'tRNA' in count_file_name:
        return pd.read_table(count_file_name,names=['id','count']) \
                .assign(type = 'tRNA') \
                .assign(name = lambda d: d['id'])
    else:
        return readGenes(count_file_name)

def readSample(count_file_path, sample_id):    
    print 'Running %s' %sample_id
    count_files = glob.glob(count_file_path + '/' + sample_id + '*counts')
    dfs = map(readDF,count_files)
    df = pd.concat(dfs, axis= 0)\
        .assign(sample_name = sample_id)
    return df

def main():
    work = os.environ['WORK']
    ref_path = os.environ['REF']
    count_file_path = '/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/mergeBam/countFiles'
    count_files = glob.glob(count_file_path + '/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))
    dfFunc = partial(readSample, count_file_path)
    dfs = Pool(12).map(dfFunc, sample_ids)
    df = pd.concat(dfs, axis=0) \
        .pipe(pd.pivot_table,index = ['id','name','type'],  
            values = 'count' , columns = ['sample_name']) \
        .fillna(0)
    tablename = count_file_path + '/combined_gene_count.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    print 'Written %s' %tablename


if __name__ == '__main__':
    main()
