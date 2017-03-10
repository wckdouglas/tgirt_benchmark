#!/bin/bash

DATAPATH=${SCRATCH}/bench_marking/data
INDEX_PATH=${REF}/human_transcriptome
INDEX=${INDEX_PATH}/transcriptome
RESULTPATH=${SCRATCH}/bench_marking/alignment_free
COUNT_PATH=${RESULTPATH}/countFiles
BAM_PATH=${RESULTPATH}/bamFiles
THREADS=24
mkdir -p ${COUNT_PATH} ${BAM_PATH}

for R1 in ${DATAPATH}/*R1_001.fastq.gz
do
	R2=${R1/R1/R2}
	SAMPLENAME=$(basename ${R1%_R1_001.fastq.gz})
	echo kallisto quant \
		-i ${INDEX} -o ${COUNT_PATH}/${SAMPLENAME} \
<<<<<<< HEAD
		--fr-stranded  --threads=${THREADS}\
		${R1} ${R2} --plainplaintext 
=======
		--bootstrap-samples=10 \
		--bias  \
		--fr-stranded --threads=${THREADS}\
		${R1} ${R2} 
>>>>>>> 6c1684cb44377ee70adabc8e4a4caf5ed6607013
done

