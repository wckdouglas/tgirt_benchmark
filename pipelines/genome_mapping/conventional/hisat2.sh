#!/bin/bash

PROJECT_PATH=${SCRATCH}/bench_marking
#PROJECT_PATH=${WORK}/cdw2854/bench_marking
PROJECT_PATH=${PROJECT_PATH}/genome_mapping
DATAPATH=${PROJECT_PATH}/Trim
RESULTPATH=${DATAPATH}/conventional
BAM_PATH=$RESULTPATH/bam_files
REFERENCE_PATH=$REF/benchmarking/GRCH38_genome
HISAT_INDEX=$REFERENCE_PATH/reference
SPLICE_FILE=$REFERENCE_PATH/splicesite.tsv
mkdir -p $BAM_PATH

for FASTQ1 in $DATAPATH/*.1.fq.gz
do
	SAMPLENAME=$(basename ${FASTQ1%.1.fq.gz})
	FASTQ2=${FASTQ1/.1./.2.}
	echo time hisat2 -p 24 -k 10 --no-mixed --no-discordant \
		--known-splicesite-infile $SPLICE_FILE  \
		-x $HISAT_INDEX -1 $FASTQ1 -2 $FASTQ2  \
		\| samtools view -bS - \
		\> $BAM_PATH/${SAMPLENAME}.bam
done



