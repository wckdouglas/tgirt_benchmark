#!/bin/bash


#tRNA.fa comes from jun

REF_PATH=/stor/work/Lambowitz/ref
TRANSCRIPTOME=$REF_PATH/human_transcriptome
ENSEMBL_GENE=ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
ENSEMBL_GTF=ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz


#Download ERCC
ERCC_annotation=https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt
curl $ERCC_annotation > $TRANSCRIPTOME/ercc_annotation.tsv
curl https://tools.thermofisher.com/content/sfs/manuals/cms_095047.txt \
	| sed 1d \
	| awk '{printf ">%s\n%s\n",$1,$5}'  \
	> $TRANSCRIPTOME/ercc.fa
cat $TRANSCRIPTOME/ercc.fa \
   | seqkit fx2tab \
   | awk '{print $1,0,length($2),$1,0,"+"}' OFS='\t' \
   > $TRANSCRIPTOME/ercc.bed

#Download transcripts and merge tRNA
curl $ENSEMBL_GENE > $TRANSCRIPTOME/transcriptome.fa.gz
zcat $TRANSCRIPTOME/transcriptome.fa.gz \
	| cat - $TRANSCRIPTOME/tRNA.fa $TRANSCRIPTOME/ercc.fa \
	> $TRANSCRIPTOME/whole_transcriptome.fa

## Download gtf
curl $ENSEMBL_GTF|gunzip > $TRANSCRIPTOME/genes.gtf


## make transcript table
OUT_FILE=$TRANSCRIPTOME/transcripts.tsv
zcat $TRANSCRIPTOME/transcriptome.fa.gz | python transcript_table_from_fa.py > $OUT_FILE
python tRNA_fai2table.py $TRANSCRIPTOME/tRNA.fa >> $OUT_FILE
awk '{print $1,$1,$1,ERCC}' OFS='\t' $TRANSCRIPTOME/ercc.bed >> $OUT_FILE
