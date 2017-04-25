#!/bin/bash

PROJECT_PATH=${SCRATCH}/bench_marking
PROJECT_PATH=${WORK}/cdw2854/bench_marking
DATAPATH=${PROJECT_PATH}/data
RESULT_PATH=${PROJECT_PATH}/genome_mapping
REFERENCE_PATH=$REF/benchmarking/genome/GRCH38_genome/
HISAT_INDEX=$REFERENCE_PATH/reference
SPLICE_FILE=$REFERENCE_PATH/splicesite.tsv
BOWTIE_INDEX=$REF/RNASeqConsortium/reference.fasta

for FQ in ${DATAPATH}/*R1_001.fastq.gz
do
	SAMPLENAME=$(basename ${FQ%_R1_001.fastq.gz})
	echo bash Hisat2_pipeline7.sh $RESULT_PATH $SAMPLENAME $HISAT_INDEX $SPLICE_FILE $BOWTIE_INDEX
done
