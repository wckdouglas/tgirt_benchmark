#!/usr/bin/env python

import pandas as pd
import glob
import os
import re
from make_customized_table import read_flag_stat, raw_read

def read_log(logfile):
    tf = []
    af = []
    for l in open(logfile):
        if 'Total fragments' in l:
            tf.append(l.split(':')[1].strip(' |\n'))
        elif 'Successfully assigned fragments' in l:
            af.append(l.split(':')[1].strip(' |\n'))
    return tf, af


def features():
    path = '/stor/home/cdw2854/tgirt_benchmark/pipelines/genome_mapping/conventional'
    log = path + '/log'
    count = '/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/Trim/conventional/counts/counts.tsv'
    header = open(count).next()
    bams = filter(lambda x: re.search('bam"$',x), header.split(' '))
    bams = pd.Series(bams).str.split('/',expand=True).values[:,-1]
    tf, af = read_log(log)
    df = pd.DataFrame({'samplename':bams,
                  'total frag': tf,
                  'Assigned frag': af}) \
        .assign(samplename = lambda d: d.samplename.str.replace('.bam"','')) 
    return df


def read_file(flag_stat_file):
    samplename = os.path.basename(flag_stat_file).replace('.flagstat','')
    fs = read_flag_stat(flag_stat_file)
    trimmed = fs['read1'][0]
    proper_map = fs['proper pair'][0]/2
    return samplename, trimmed, proper_map


project_path = '/stor/work/Lambowitz/cdw2854/bench_marking'
base_path = project_path + '/genome_mapping/Trim'
flagstat_path = base_path + '/conventional/bam_files'
samples = glob.glob(flagstat_path + '/*flagstat')
df = pd.DataFrame(map(read_file, samples), 
                  columns = ['samplename','trimmed','proper mapped'])
feature_df = features()
df = df.merge(feature_df, on='samplename') \
        .merge(raw_read(), on='samplename')\
        .sort_values('samplename') \
        .pipe(lambda d: d[['samplename','raw count','trimmed','proper mapped',
                           'Assigned frag']]) 
df['Assigned frag'] = df['Assigned frag'].str.split(' ', expand=True)[0]
df['proper mapped'] = df['proper mapped'] * 2
df.columns = ['Sample name','Raw pairs','Trimmed pairs','Properly mapped fragments','Feature fragments']
df.columns = df.columns.str.replace('_',' ')
tablename = project_path + '/tables/conventional_table.csv'
df.to_csv(tablename, index=False)
print 'Written %s' %tablename




