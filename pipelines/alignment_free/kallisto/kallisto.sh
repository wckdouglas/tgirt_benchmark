#!/bin/bash

DATAPATH=${SCRATCH}/bench_marking/genome_mapping/Trim
INDEX_PATH=${REF}/benchmarking/human_transcriptome
INDEX=${INDEX_PATH}/transcriptome_kallisto/transcriptome_kallisto
RESULTPATH=${SCRATCH}/bench_marking/alignment_free
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

