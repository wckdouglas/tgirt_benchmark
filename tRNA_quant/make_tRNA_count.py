#!/usr/bin/env python

import pysam
import glob
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import os

def count_tRNA(bam_file):
    samplename = os.path.basename(bam_file).split('.')[0]
    print 'Counting %s' %samplename
    tRNA_count = defaultdict(int)
    with pysam.Samfile(bam_file) as bam:
        for aln in bam:
            if aln.flag == 0:
                tRNA = aln.reference_name
                tRNA_count[tRNA] += 1

    print 'Finished counting %s' %samplename
    return pd.DataFrame({'tRNA': tRNA_count.keys(),
                  'read_count':tRNA_count.values(),
                  'samplename':samplename})
                

def main():
    project_path ='/stor/work/Lambowitz/cdw2854/bench_marking/tRNA_seq/'
    bam_path = project_path + '/bam_files'
    count_path = project_path + '/tRNA_counts'

    if not os.path.isdir(count_path):
        os.mkdir(count_path)

    bam_files = glob.glob(bam_path+ '/*bam')
    p = Pool(8)
    df = p.map(count_tRNA, bam_files)
    p.close()
    p.join()

    df = pd.concat(df, axis=0)
    tablename = count_path + '/count_table.tsv'
    df.to_csv(tablename, sep='\t',index=False)
    print 'Written', tablename


if __name__ == '__main__':
    main()
