#!/bin/bash

REF_PATH=$REF/benchmarking/human_transcriptome
TRNA_FA=$REF_PATH/tRNA.fa
TRANSCRIPTS_BED=$REF_PATH/genes.bed
TRANSCRIPT_LENGTH=$REF_PATH/genes.length

cat $TRANSCRIPTS_BED \
	| awk '{print $NF, $3-$2}' OFS='\t'  \
	| sed 1i"id\tgene_length"  \
	> $TRANSCRIPT_LENGTH

samtools faidx $TRNA_FA
cat ${TRNA_FA}.fai \
	| cut -f1,2 \
	>> $TRANSCRIPT_LENGTH	
