#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import os
import re
from builtins import map, range

def trim():
    seqfile = '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking/genome_mapping/tgirt_map/Trim/seq_count.txt'
    df = pd.read_table(seqfile, sep=':', names=['samplename','trim_count']) \
        .assign(samplename = lambda d: d.samplename.str.replace('.1.fq.gz','').map(os.path.basename)) 
    return df

def raw_read():
    seq_count = '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking/data/seq_count.txt'
    seq_df = pd.read_table(seq_count, sep=':',names=['samplename','raw count'])\
            .assign(samplename = lambda d: d.samplename.str.replace('.fastq.gz','').map(os.path.basename))
    return seq_df


def get_number(line):
    return [int(line.split(' ')[0])]

def get_line(lines, keyword):
    return_line = ''
    for line in lines:
        if keyword in line:
            return_line += line
    return return_line

def read_flag_stat(filename):
    info = open(filename,'r').readlines()
    stat_df = {}
    stat_df['read1'] = get_number(get_line(info, 'read1'))
    stat_df['mapped'] = get_number(get_line(info, 'mapped'))
    stat_df['supplementary'] = get_number(get_line(info, 'supplementary'))
    stat_df['proper pair'] = get_number(get_line(info, 'properly paired'))
    stat_df['secondary'] = get_number(get_line(info, 'secondary'))
    return stat_df


def read_sample(sample_folder):
    samplename = os.path.basename(sample_folder)
    hisat = read_flag_stat(sample_folder + '/Hisat/hisat.flagstat')
    trimmed_non_trRNA = hisat['read1'][0]
    hisat_mapped = hisat['proper pair'][0]/2

    hisat_unique_mapped = read_flag_stat(sample_folder+'/Hisat/hisat.unique.flagstat')['read1'][0]
    
    bowtie = read_flag_stat(sample_folder+'/Bowtie/bowtie2.flagstat')
    bowtie_mapped = bowtie['proper pair'][0]/2 
    bowtie_unique_mapped = read_flag_stat(sample_folder + '/Bowtie/bowtie.unique.flagstat')['read1'][0]

    premap = read_flag_stat(sample_folder + '/rRNA_tRNA_premap/tRNA_rRNA.flagstat')
    premap_reads = premap['mapped'][0]/2
    return (samplename, premap_reads, trimmed_non_trRNA,  
            hisat_mapped, hisat_unique_mapped, 
            bowtie_mapped, bowtie_unique_mapped)

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'
    base_path = project_path + '/genome_mapping/tgirt_map'
    samples = glob.glob(base_path + '/Sample-*')
    samples = filter(lambda x: re.search('-[ABCD]-[123]',x), samples)
    df = pd.DataFrame(list(map(read_sample,samples)), 
                      columns = [ 'samplename','premap_trRNA','trimmed_non_trRNA',
                                'hisat_mapped','hisat_uniq','bowtie_mapped','bowtie_uniq'])
    df = df.merge(trim()).merge(raw_read()) \
        .sort_values('samplename')\
        .pipe(lambda d: d[['samplename','raw count','trim_count','premap_trRNA','trimmed_non_trRNA','hisat_mapped',
             'hisat_uniq','bowtie_mapped','bowtie_uniq']]) \
        .assign(samplename = lambda d: d.samplename.str.replace('_R1_001$',''))
    df.columns = ['Sample name','Raw pairs','Trimmed pairs','tRNA/rRNA pairs','non tRNA/rRNA pairs','HISAT mapped pairs',
                  'HISAT uniquely mapped','BOWTIE2 mapped pairs','BOWTIE2 uniquely mapped']
    tablename = project_path + '/tables/customized_table.csv'
    df.to_csv(tablename, index=False)
    print('Written %s' %tablename)

if __name__ == '__main__':
    main()
