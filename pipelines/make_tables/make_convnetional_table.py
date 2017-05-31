#!/usr/bin/env python

import pandas as pd
import glob
import os
from make_customized_table import read_flag_stat



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
print flagstat_path
df = pd.DataFrame(map(read_file, samples), 
                  columns = ['samplename','trimmed','proper mapped'])
tablename = project_path + '/tables/conventional_table.tsv'
df.to_csv(tablename, sep='\t', index=False)
print 'Written %s' %tablename




