#!/bin/bash

DATAPATH=${SCRATCH}/bench_marking/data
INDEX_PATH=${REF}/human_transcriptome
INDEX=${INDEX_PATH}/transcript_salmon
RESULTPATH=${SCRATCH}/bench_marking/alignment_free/salmon
COUNT_PATH=${RESULTPATH}/countFiles
BAM_PATH=${RESULTPATH}/bamFiles
THREADS=24
mkdir -p ${COUNT_PATH} ${BAM_PATH}

for R1 in ${DATAPATH}/*R1_001.fastq.gz
do
	R2=${R1/R1/R2}
	SAMPLENAME=$(basename ${R1%_R1_001.fastq.gz})
	echo salmon quant \
		--seqBias --gcBias \
		--index $INDEX --libType ISF \
		--writeMappings  \
		--threads=${THREADS} \
		--auxDir $RESULTPATH/${SAMPLENAME} \
		--numBootstraps 100 \
		--mates1 ${R1} --mates2 ${R2} \
		--output $RESULTPATH/${SAMPLENAME} \
		\| samtools view -b@ $THREADS \
		\> ${RESULTPATH}/${SAMPLENAME}.bam
done

