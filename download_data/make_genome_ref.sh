#!/bin/bash

REF_PATH=$REF/RNASeqConsortium
TRANSCRIPTOME_PATH=$REF/human_transcriptome
HUMAN_REF_URL=ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
curl $HUMAN_REF_URL > $REF_PATH/reference.fa.gz
zcat $REF_PATH/reference.fa.gz \
	| cat - $TRANSCRIPTOME_PATH/ercc.fa $TRANSCRIPTOME_PATH/rRNA.fa \
	> reference.fa 
cd $REF_PATH
hisat2-build reference.fa reference.fasta
