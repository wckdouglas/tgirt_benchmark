#!/usr/bin/env python

from collections import defaultdict

raw_data_path = '/stor/scratch/Lambowitz/cdw2854/bench_marking/data'

sample_dict = defaultdict(int)
with open('SraRunTable.txt','r') as table:
    for i, line in enumerate(table):
        if i!=0:
            fields = line.split('\t')
            samplename = fields[2].replace(' ','-')
            run_number = fields[6]
            sample_dict[samplename] += 1
            current = sample_dict[samplename]
            sample_id =  samplename+'-'+str(current) 

            #render commands
            download_command = 'fastq-dump --split-3 --gzip --outdir %s %s ' %(raw_data_path, run_number)
            rename_read1_command = 'mv %s/%s_1.fastq.gz %s/%s_R1_001.fastq.gz' \
                    %(raw_data_path,run_number, raw_data_path, sample_id)
            rename_read2_command = 'mv %s/%s_2.fastq.gz %s/%s_R2_001.fastq.gz' \
                    %(raw_data_path,run_number, raw_data_path, sample_id)
            print ';'.join([download_command, rename_read1_command, rename_read2_command])
