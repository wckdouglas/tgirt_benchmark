#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking
DATAPATH=$PROJECT_PATH/genome_mapping/tgirt_map/Trim
REF_PATH=$REF
INDEX_PATH=${REF_PATH}/benchmarking/human_transcriptome
INDEX=${INDEX_PATH}/transcript_salmon
RESULTPATH=$PROJECT_PATH/alignment_free/salmon
THREADS=24

for R1 in ${DATAPATH}/*.1.fq.gz
do
	for KMER in '' '_11' '_15' '_21'
	do
        BAM_PATH=${RESULTPATH}${KMER}/bamFiles
        mkdir -p ${BAM_PATH}
		R2=${R1/.1./.2.}
		SAMPLENAME=$(basename ${R1%.1.fq.gz})
        OUT_PATH=${RESULTPATH}${KMER}/${SAMPLENAME}
		echo time salmon quant \
			--seqBias --gcBias \
			--index ${INDEX}${KMER} --libType ISF \
			--writeMappings  \
			--threads=${THREADS} \
			--auxDir aux \
			--numBootstraps 100 \
			--mates1 ${R1} --mates2 ${R2} \
			--output $OUT_PATH \
			\| samtools view -b@ $THREADS \
			\> ${BAM_PATH}/${SAMPLENAME}.bam
	done
done

