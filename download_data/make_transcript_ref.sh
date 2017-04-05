#!/bin/bash


#tRNA.fa comes from jun

REF_PATH=/stor/work/Lambowitz/ref
TRANSCRIPTOME=$REF_PATH/human_transcriptome
ENSEMBL_TRANSCRIPT=ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
ENSEMBL_NON_CODING=ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
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
   | awk '{print $1,0,length($2),$1,0,"+","ERCC",$1}' OFS='\t' \
   > $TRANSCRIPTOME/ercc.bed

#Download ERCC
python get_rRNA.py $TRANSCRIPTOME/rRNA
echo 'gi|23898|emb|X12811.1|  274     394     5S_rRNA 0       +       5S_rRNA 5S_rRNA
gi|555853|gb|U13369.1|HSU13369  3657    5527    18S_rRNA        0       +       18S_rRNA        18S_rRNA
gi|555853|gb|U13369.1|HSU13369  6623    6779    5.8S_rRNA       0       +       5.8S_rRNA       5.8S_rRNA
gi|555853|gb|U13369.1|HSU13369  7935    12969   28S_rRNA        0       +       28S_rRNA        28S_rRNA' > $TRANSCRIPTOME/rRNA.bed

#Download transcripts and merge tRNA
curl $ENSEMBL_TRANSCRIPT > $TRANSCRIPTOME/ensembl_cDNA.fa.gz
curl $ENSEMBL_NON_CODING > $TRANSCRIPTOME/ensembl_ncrna.fa.gz
zcat $TRANSCRIPTOME/transcriptome.fa.gz \
		$TRANSCRIPTOME/ensembl_ncrna.fa.gz \
	| python correct_transcriptome_id.py \
	| tee $TRANSCRIPTOME/ensembl_transcripts.fa \
	| cat - $TRANSCRIPTOME/tRNA.fa $TRANSCRIPTOME/rRNA.fa $TRANSCRIPTOME/ercc.fa \
	> $TRANSCRIPTOME/whole_transcriptome.fa

## make transcript table
OUT_FILE=$TRANSCRIPTOME/transcripts.tsv
cat $TRANSCRIPTOME/ensembl_transcripts.fa \
	| python transcript_table_from_fa.py > $OUT_FILE
python tRNA_fai2table.py $TRANSCRIPTOME/tRNA.fa >> $OUT_FILE
awk '{print $1,$1,$1,"ERCC"}' OFS='\t' $TRANSCRIPTOME/ercc.bed >> $OUT_FILE
awk -F'\t' '{print $1,$4,$1,"rRNA"}' OFS='\t' $TRANSCRIPTOME/rRNA.bed >> $OUT_FILE
