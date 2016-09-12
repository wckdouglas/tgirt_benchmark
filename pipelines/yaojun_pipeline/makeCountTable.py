#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import glob

def readGenes(count_file_name):
    return pd.read_table(count_file_name,names = [''])

def readSample(count_file_path, sample_id):    
    count_files = glob.glob(count_file_path + '/' + sample_id + '*counts')


def readRibosomal(ref_path):
    bed_path = ref_path + '/GRCh38/Bed_for_counts_only'
    rDNA = pd.read_table(bed_path + '/rDNA.bed', header=None)[[0,3,6]] 
    rDNA.columns = ['id','type','name']
    return rDNA



def main():
    work = os.environ['WORK']
    ref_path = os.environ['REF']
    count_file_path = '/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/mergeBam/counts'
    count_files = glob.glob(count_file_path + '/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))

    rDNA_genes = readRibosomal(ref_path)





if __name__ == '__main__':
    main()
