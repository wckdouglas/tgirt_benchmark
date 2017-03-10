#!/bin/bash

PROJECT_PATH=${SCRATCH}/bench_marking
DATAPATH=${PROJECT_PATH}/data
RESULT_PATH=${PROJECT_PATH}/genome_mapping
REF_PATH=${REF}
HISAT_INDEX=$REF/GRCh38/Plasma_ref/hisat2 
SPLICE_FILE=$REF/GRCh38/Plasma_ref/splicesites.txt
BOWTIE_INDEX=$REF/GRCh38/Plasma_ref/plasma_ref

for FQ in ${DATAPATH}/*R1_001.fastq.gz
do
	SAMPLENAME=$(basename ${FQ%_R1_001.fastq.gz})
	echo bash Hisat2_pipeline5.sh $RESULT_PATH $SAMPLENAME $HISAT_INDEX $SPLICE_FILE $BOWTIE_INDEX
done
