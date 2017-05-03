#!/usr/bin/env python

import os
import glob
import re

work_path = os.environ['WORK']
scratch_path = os.environ['SCRATCH']
gene_file = os.environ['REF'] + '/benchmarking/human_transcriptome/genes.SAF'
project_path = work_path + '/cdw2854/bench_marking/genome_mapping'
map_path = project_path + '/Trim/conventional/bam_files'
count_path = project_path + '/Trim/conventional/counts'
if not os.path.isdir(count_path):
    os.mkdir(count_path)
bam_files = glob.glob(map_path + '/*bam') 
samplename = map(lambda x: os.path.basename(x).replace('.bam',''), bam_files)

feature_count_command = 'featureCounts -a %s ' %(gene_file) +\
                        '%s ' %(' '.join(bam_files)) + \
                        '-F SAF ' +\
                        '-O ' +\
                        '-s 1 ' +\
                        '-M ' +\
                        '-T 24 '+\
                        '--largestOverlap ' +\
                        '--minOverlap 10 '+ \
                        '--primary '+\
                        '-p ' +\
                        '-P ' +\
                        '-d 10 '+\
                        '-D 10000 '+\
                        '-B ' +\
                        '-S fr ' +\
                        '-C '+\
                        '--donotsort '+\
                        '-o /dev/stdout >  %s/counts.tsv ' %(count_path) 
print 'cd %s' %count_path
print feature_count_command
