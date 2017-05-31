#!/bin/bash

DATAPATH=${SCRATCH}/bench_marking/genome_mapping/Trim
INDEX_PATH=${REF}/benchmarking/human_transcriptome
INDEX=${INDEX_PATH}/transcript_salmon
RESULTPATH=${SCRATCH}/bench_marking/alignment_free/salmon
COUNT_PATH=${RESULTPATH}/countFiles
BAM_PATH=${RESULTPATH}/bamFiles
THREADS=24
mkdir -p ${COUNT_PATH} ${BAM_PATH}

for R1 in ${DATAPATH}/*.1.fq.gz
do
	R2=${R1/.1./.2.}
	SAMPLENAME=$(basename ${R1%.1.fq.gz})
	echo time salmon quant \
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

