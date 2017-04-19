#!/bin/bash

TRNA_FA=/stor/work/Lambowitz/ref/human_transcriptome/tRNA.fa.fai
TRANSCRIPTS_BED=/stor/work/Lambowitz/ref/RNASeqConsortium/genes.bed
TRANSCRIPT_LENGTH=/stor/work/Lambowitz/ref/human_transcriptome/genes.length

cat $TRANSCRIPTS_BED \
	| awk '{print $NF, $3-$2}' OFS='\t'  \
	| sed 1i"id\tgene_length"  \
	> $TRANSCRIPT_LENGTH

cat $TRNA_FA \
	| cut -f1,2 \
	>> $TRANSCRIPT_LENGTH	
