#!/usr/bin/env python

import os
import glob
import re

work_path = os.environ['WORK']
scratch_path = os.environ['SCRATCH']
gene_file = os.environ['REF'] + '/human_transcriptome/genes.gtf'
project_path = work_path + '/cdw2854/bench_marking/genome_mapping/stringtie'
project_path = scratch_path + '/bench_marking/genome_mapping/stringtie'
bam_path = project_path + '/bam_files' 
count_path = project_path + '/counts'
if not os.path.isdir(count_path):
    os.mkdir(count_path)
bam_files = glob.glob(bam_path + '/*.bam') 

for bam_file in bam_files:
    samplename = os.path.basename(bam_file).split('.')[0]
    result_path = 
    command = 'stringtie -p 12 '+\
            '-G %s ' %(gene_file)+\
            '-o %s/%s.gtf ' %(count_path, samplename) +\
            '-A %s/%s.gene_counts ' %(count_path, samplename) +\
            '%s ' %(bam_file)
    print command
