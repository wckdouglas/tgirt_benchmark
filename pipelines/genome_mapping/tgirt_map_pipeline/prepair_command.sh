#!/bin/bash

PROJECT_PATH=${SCRATCH}/bench_marking
DATAPATH=${PROJECT_PATH}/data
RESULT_PATH=${PROJECT_PATH}/genome_mapping/tgirt_map
REF_PATH=$SCRATCH/ref/benchmarking_new
HUMAN_INDEX=${REF_PATH}/benchmarking/GRCH38_genome/reference  # GRCh38 hisat2 index 
BED_PATH=${REF_PATH}/benchmarking/human_transcriptome  #bed file GRCh38
tRNA_INDEX=${REF_PATH}/benchmarking/human_transcriptome/tRNA #tRNA with mtTRNA and cytosolic tRNA
rRNA_INDEX=${REF_PATH}/benchmarking/human_transcriptome/rRNA #tRNA with mtTRNA and cytosolic tRNA
trRNA_INDEX=${REF_PATH}/benchmarking/human_transcriptome/tRNA_rRNA #tRNA with mtTRNA and cytosolic tRNA
SPLICE_FILE=${REF_PATH}/benchmarking/GRCH38_genome/splicesite.txt
THREADS=24
mkdir -p $RESULT_PATH

for FQ1 in ${DATAPATH}/*R1_001.fastq.gz
do
	FQ2=${FQ1/R1/R2}
    echo tgirt_count.py -1 ${FQ1} -2 ${FQ2}\
			-o ${RESULT_PATH} \
			-x ${HUMAN_INDEX} \
			-y ${HUMAN_INDEX} \
			-b ${BED_PATH} \
		    -s ${SPLICE_FILE} \
			-t ${tRNA_INDEX} \
			-r ${rRNA_INDEX} \
	    	-p ${THREADS} \
			-e ${trRNA_INDEX}
done
