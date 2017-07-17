#!/bin/bash

tRNA_PATH=/stor/work/Lambowitz/cdw2854/bench_marking/tRNA_seq
DATA_PATH=$tRNA_PATH/data
DATA_PATH=/stor/work/Lambowitz/Data/NGS/Matt
BAM_PATH=$tRNA_PATH/bam_files
REF_PATH=$REF/benchmarking/human_transcriptome/tRNA

for FQ in $DATA_PATH/*fastq.gz
do
	SAMPLENAME=$(basename ${FQ%.fastq.gz})
	echo bowtie2 --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 \
			-k 10 -x $REF_PATH -U $FQ  \
		\| samtools view -b \
		\> $BAM_PATH/${SAMPLENAME}.bam
done
