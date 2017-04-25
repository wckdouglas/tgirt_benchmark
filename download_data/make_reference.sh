#!/bin/bash

REF_PATH=${REF}/benchmarking
TRANSCRIPTOME=$REF_PATH/human_transcriptome
GENOME_PATH=$REF_PATH/GRCH38_genome
mkdir -p $TRANSCRIPTOME $GENOME_PATH

#URLs for reference
ENSEMBL_TRANSCRIPT=ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
ENSEMBL_NON_CODING=ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
ENSEMBL_GTF=ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.chr_patch_hapl_scaff.gtf.gz
HUMAN_REF_URL=ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ERCC_annotation=https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt
GTRNA=http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz

### Download genome ref
#curl $HUMAN_REF_URL > $GENOME_PATH/reference.fa.gz

#Download rRNA
python get_rRNA_fa.py > $TRANSCRIPTOME/rRNA.fa
echo 'gi|23898|emb|X12811.1|  274     394     5S_rRNA 0       +       5S_rRNA 5S_rRNA
gi|555853|gb|U13369.1|HSU13369  3657    5527    18S_rRNA        0       +       18S_rRNA        18S_rRNA
gi|555853|gb|U13369.1|HSU13369  6623    6779    5.8S_rRNA       0       +       5.8S_rRNA       5.8S_rRNA
gi|555853|gb|U13369.1|HSU13369  7935    12969   28S_rRNA        0       +       28S_rRNA        28S_rRNA' \
	| awk '{print $1,$2,$3,$4,$5,$6,"rDNA",$8}' OFS='\t' \
	> $TRANSCRIPTOME/rRNA.bed
echo 'Made rRNA'

#Download ERCC
curl $ERCC_annotation | python clean_ercc.py  > $TRANSCRIPTOME/ercc_annotation.tsv
curl https://tools.thermofisher.com/content/sfs/manuals/cms_095047.txt \
	| sed 1d \
	| awk '{printf ">%s\n%s\n",$1,$5}'  \
	> $TRANSCRIPTOME/ercc.fa
cat $TRANSCRIPTOME/ercc.fa \
   | seqkit fx2tab \
   | awk '{print $1,0,length($2),$1,0,"+","ERCC",$1}' OFS='\t' \
   > $TRANSCRIPTOME/ercc.bed
echo 'Made ERCC'

####MAKE hisat2 index
zcat $GENOME_PATH/reference.fa.gz \
    | cat - $TRANSCRIPTOME/ercc.fa $TRANSCRIPTOME/rRNA.fa \
> $GENOME_PATH/reference.fa
samtools faidx $GENOME_PATH/reference.fa
#hisat2-build $GENOME_PATH/reference.fa $GENOME_PATH/reference
#bowtie2-build $GENOME_PATH/reference.fa $GENOME_PATH/reference
hisat2_extract_splice_sites.py $GENOME_PATH/genes.gtf > $GENOME_PATH/splicesite.tsv  
echo 'Made genome'

#download gene annotation
Rscript get_gene_bed.R $TRANSCRIPTOME/genes.bed
cat $TRANSCRIPTOME/rRNA.bed >> $TRANSCRIPTOME/genes.bed
cat $TRANSCRIPTOME/ercc.bed >> $TRANSCRIPTOME/genes.bed
echo 'Made genes.bed'

#download tRNA
tRNA_PATH=$TRANSCRIPTOME/tRNA
mkdir -p $tRNA_PATH
curl -o $tRNA_PATH/hg38-tRNAs.tar.gz $GTRNA  
tar zxvf $tRNA_PATH/hg38-tRNAs.tar.gz --directory $tRNA_PATH
python scrape_tRNA_name.py $tRNA_PATH
python make_tRNA_fasta.py $tRNA_PATH > $tRNA_PATH/nucleo_tRNA.fa
cat $tRNA_PATH/hg38_tRNA.info | sed 1d > $TRANSCRIPTOME/tRNA.bed 
cat $TRANSCRIPTOME/tRNA.bed |cut -f1-8 >> $TRANSCRIPTOME/genes.bed
cat $TRANSCRIPTOME/genes.bed \
	| grep 'Mt_tRNA' \
	| bedtools getfasta  -fi $GENOME_PATH/reference.fa -bed - -s -name -tab \
	| tr ':' '\t' \
	| awk '{printf ">%s\n%s\n",$1,$NF}' \
	>> $tRNA_PATH/mt_tRNA.fa
echo 'Finished making tRNA'
cat $tRNA_PATH/mt_tRNA.fa $tRNA_PATH/nucleo_tRNA.fa > $TRANSCRIPTOME/tRNA.fa


#Download transcripts and merge tRNA, ercc, rDNA
curl $ENSEMBL_TRANSCRIPT > $TRANSCRIPTOME/ensembl_cDNA.fa.gz
curl $ENSEMBL_NON_CODING > $TRANSCRIPTOME/ensembl_ncrna.fa.gz
bedtools getfasta -s -fi $TRANSCRIPTOME/rRNA.fa -bed $TRANSCRIPTOME/rRNA.bed -name > $TRANSCRIPTOME/rDNA.fa
zcat $TRANSCRIPTOME/ensembl_cDNA.fa.gz \
		$TRANSCRIPTOME/ensembl_ncrna.fa.gz \
	| python correct_transcriptome_id.py \
	| tee $TRANSCRIPTOME/ensembl_transcripts.fa \
	| cat - $tRNA_PATH/nucleo_tRNA.fa $TRANSCRIPTOME/rDNA.fa $TRANSCRIPTOME/ercc.fa \
	> $TRANSCRIPTOME/whole_transcriptome.fa
samtools faidx $TRANSCRIPTOME/whole_transcriptome.fa
echo 'Made transcriptome fasta'

## make transcript table
OUT_FILE=$TRANSCRIPTOME/transcripts.tsv
cat $TRANSCRIPTOME/ensembl_transcripts.fa \
	| python transcript_table_from_fa.py > $OUT_FILE
awk '{print $1,$1,$1,"ERCC"}' OFS='\t' $TRANSCRIPTOME/ercc.bed >> $OUT_FILE
awk -F'\t' '{print $1,$4,$1,"rRNA"}' OFS='\t' $TRANSCRIPTOME/rRNA.bed >> $OUT_FILE
cat $tRNA_PATH/hg38_tRNA.info \
	| awk '{print $8, $4, $8, $7}' OFS='\t'\
	| sed 1d    \
	>> $OUT_FILE
echo 'Made transcript table'

bowtie2-build $TRANSCRIPTOME/tRNA.fa $TRANSCRIPTOME/tRNA
bowtie2-build $TRANSCRIPTOME/rRNA.fa $TRANSCRIPTOME/rRNA

## Download GTF and append tRNA, rRNA, ERCC bed record
GENES_GTF=$GENOME_PATH/genes.gtf
curl $ENSEMBL_GTF  \
	| zcat \
	| grep -v 'gene_biotype "TEC"' \
	> $GENES_GTF
python bed_to_gtf.py $TRANSCRIPTOME/rRNA.bed >> $GENES_GTF
python bed_to_gtf.py $TRANSCRIPTOME/ercc.bed >> $GENES_GTF
python bed_to_gtf.py $TRANSCRIPTOME/tRNA.bed >> $GENES_GTF
echo 'Made gtf'
	
#for genome count 
SAF_FILE=$TRANSCRIPTOME/genes.SAF
echo 'GeneID\tchr\tstart\tend\tstrand' > $SAF_FILE
awk '{print $NF,$1,$2,$3,$6}' OFS='\t' $TRANSCRIPTOME/genes.bed >> $SAF_FILE
python split_bed_for_count.py $TRANSCRIPTOME $TRANSCRIPTOME/transcripts.tsv
echo 'Made bed for count and SAF file'

# get union gene set
python calibrate_gene_set.py $TRANSCRIPTOME


#make tRNA and rRNA fasta
cat $TRANSCRIPTOME/genes.bed \
	| awk '$7=="rRNA"' \
	| bedtools getfasta -fi $GENOME_PATH/reference.fa -bed -  -s \
	| fastx_collapser \
	| sed 's/>/>rRNA_/g' \
	| cat - $TRANSCRIPTOME/tRNA.fa $TRANSCRIPTOME/rRNA.fa \
	> $TRANSCRIPTOME/tRNA_rRNA.fa
bowtie2-build $TRANSCRIPTOME/tRNA_rRNA.fa $TRANSCRIPTOME/tRNA_rRNA

