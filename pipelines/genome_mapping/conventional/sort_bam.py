#!/usr/bin/env python

import os
import glob
import re

work_path = os.environ['WORK']
scratch_path = os.environ['SCRATCH']
gene_file = os.environ['REF'] + '/human_transcriptome/genes.gtf'
project_path = work_path + '/cdw2854/bench_marking/genome_mapping'
project_path = scratch_path + '/bench_marking/genome_mapping'
bam_paths = project_path + '/pipeline7'
count_path = project_path + '/stringtie/counts'
out_bam_path = project_path + '/stringtie/bam_files'
if not os.path.isdir(count_path):
    os.makedirs(count_path)
if not os.path.isdir(out_bam_path):
    os.makedirs(out_bam_path)
sample_folders = glob.glob(bam_paths + '/*') 
sample_folders = filter(lambda x: re.search('-[0-9]$',x), sample_folders)

for sample_folder in sample_folders:
    bam_file = sample_folder +'/Hisat/hisat.bam'
    samplename = os.path.basename(sample_folder)
    command = 'sambamba sort --nthreads=16 '+\
            '--tmpdir=%s/%s ' %(out_bam_path, samplename)+\
            '--out=%s/%s.hisat.sorted.bam ' %(out_bam_path, samplename) +\
            '--show-progress '+\
            '%s ' %(bam_file) 
    print command
