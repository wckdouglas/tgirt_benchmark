#!/bin/bash

PROJECT_PATH=${SCRATCH}/bench_marking
DATAPATH=${PROJECT_PATH}/data
RESULT_PATH=${PROJECT_PATH}/genome_mapping
REF_PATH=${REF}
HUMAN_INDEX=${REF}/RNASeqConsortium/reference.fasta # GRCh38 hisat2 index 
BED_PATH=${REF}/GRCh38/Bed_for_counts_only #bed file GRCh38
tRNA_INDEX=${REF}/GRCh38/tRNA/temp/tRNA #tRNA with mtTRNA and cytosolic tRNA
rRNA_INDEX=${REF}/GRCh38/Plasma_ref/rDNA #tRNA with mtTRNA and cytosolic tRNA
SPLICE_FILE=${REF}/RNASeqConsortium/splicesite.txt
THREADS=16
ADAPTOR=adaptors.fa
STRAND=forward

for FQ in ${DATAPATH}/*R1_001.fastq.gz
do
    echo time python -u tgirtPipelinePaired.py --fastq=${FQ} \
			--outdir=${RESULT_PATH} \
			--humanIndex=${HUMAN_INDEX} \
		    --splicesite=${SPLICE_FILE} \
			--tRNAindex=${tRNA_INDEX} \
	    	--threads=${THREADS} \
            --adaptors=${ADAPTOR} \
			--strand=${STRAND} \
			--bedpath=${BED_PATH} \
			--rRNAindex=${rRNA_INDEX}
done
