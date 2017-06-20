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
	for KMER in '' '_11' '_15' '_21'
	do
		R2=${R1/.1./.2.}
		SAMPLENAME=$(basename ${R1%.1.fq.gz})
		echo time salmon quant \
			--seqBias --gcBias \
			--index ${INDEX}${KMER} --libType ISF \
			--writeMappings  \
			--threads=${THREADS} \
			--auxDir ${RESULTPATH}${KMER}/${SAMPLENAME} \
			--numBootstraps 100 \
			--mates1 ${R1} --mates2 ${R2} \
			--output ${RESULTPATH}${KMER}/${SAMPLENAME} \
			\| samtools view -b@ $THREADS \
			\> ${RESULTPATH}${KMER}/${SAMPLENAME}.bam
	done
done

