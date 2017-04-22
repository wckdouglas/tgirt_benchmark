#!/bin/bash

REF_PATH=$REF/benchmarking/GRCH38_genome
mkdir -p $REF_PATH
TRANSCRIPTOME_PATH=$REF/benchmarking/human_transcriptome

#download fasta
HUMAN_REF_URL=ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
#curl $HUMAN_REF_URL > $REF_PATH/reference.fa.gz
#gunzip $REF_PATH/reference.fa.gz
cat $REF_PATH/reference.fa \
	| cat - $TRANSCRIPTOME_PATH/ercc.fa $TRANSCRIPTOME_PATH/rRNA.fa \
> reference.fa 
cd $REF_PATH
hisat2-build reference.fa reference.fasta
