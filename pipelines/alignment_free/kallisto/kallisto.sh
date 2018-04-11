#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking
DATAPATH=${PROJECT_PATH}/genome_mapping/tgirt_map/Trim
INDEX_PATH=${REF}/benchmarking/human_transcriptome
INDEX=${INDEX_PATH}/transcriptome_kallisto/transcriptome_kallisto
RESULTPATH=${PROJECT_PATH}/alignment_free
COUNT_PATH=${RESULTPATH}/kallisto
BAM_PATH=${RESULTPATH}/bamFiles
THREADS=24
mkdir -p ${COUNT_PATH} ${BAM_PATH}

for R1 in ${DATAPATH}/*.1.fq.gz
do
	R2=${R1/.1./.2.}
	SAMPLENAME=$(basename ${R1%.1.fq.gz})
	echo time kallisto quant \
		-i ${INDEX} -o ${COUNT_PATH}/${SAMPLENAME} \
		--fr-stranded  --threads=${THREADS}\
		${R1} ${R2}  
done

