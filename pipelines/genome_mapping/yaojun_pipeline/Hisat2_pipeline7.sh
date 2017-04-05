#!bin/bash

# Pipeline for PE reads

# Need JA number, and part of the sample name from the RAW fastq file as input parameters
# Such as JA15599, and YQ1
# $1 is the JA number start with JA, and $2 is the sample name, better to choose a short but not an ambiguous one
# $3 is the reference genome path, the index prefix
# $4 is the known splice site file, in txt format
# $5 is the Bowtie2 index

# 1. check folder structure, if not exit, make one
# folder structure:
# RAW fastq is stored under $HOME/NGS/Data/$1, which should exist before running this script
# The folder will be checked and created while not existing is $HOME/NGS/Work/$1
# $HOME/NGS/Work/$1-----Trim
# 		  |-----$2 (your sample name)-----Hisat (Pass1)
# 		  |			    |-----Bowtie (Pass2)
# 		  |			    |-----Combined
# 		  |		            |-----tRNA
#                 |                         |-----rRNA
# 		  |-----Counts-----RAW
# 		             |-----Simple
# 		             |-----tRNA_RAW
# 		             |-----tRNA_anti
Workfolder=$1
Sample=$2
Ref=$3
Splicesite=$4
Bowtie=$5

RAWfolder=$(dirname $Workfolder)/data
Trimfolder=$Workfolder/Trim
Samplefolder=$Workfolder/$Sample
Countsfolder=$Workfolder/Counts

if [ ! -d "$RAWfolder" ] ; then echo "This folder doesn't exist."; exit 1; fi
if [ ! -d "$Workfolder" ] ; then
	mkdir -p "$Trimfolder"
	mkdir -p "$Samplefolder/Hisat"
	mkdir "$Samplefolder/Bowtie"
	mkdir "$Samplefolder/Combined"
	mkdir "$Samplefolder/tRNA"
	mkdir "$Samplefolder/rRNA"
	mkdir -p "$Countsfolder/RAW"
	mkdir "$Countsfolder/Simple"
	mkdir "$Countsfolder/tRNA_RAW"
	mkdir "$Countsfolder/tRNA_anti"
else
	if [ ! -d "$Samplefolder" ] ; then
		mkdir -p "$Samplefolder/Hisat"
		mkdir "$Samplefolder/Bowtie"
		mkdir "$Samplefolder/Combined"
		mkdir "$Samplefolder/tRNA"
		mkdir "$Samplefolder/rRNA"
	else
		if [ ! -d "$Samplefolder/Hisat" ] ; then
			mkdir "$Samplefolder/Hisat"
		fi
		if [ ! -d "$Samplefolder/Bowtie" ] ; then
			mkdir "$Samplefolder/Bowtie"
		fi
		if [ ! -d "$Samplefolder/Combined" ] ; then
			mkdir "$Samplefolder/Combined"
		fi
		if [ ! -d "$Samplefolder/tRNA" ] ; then
			mkdir "$Samplefolder/tRNA"
		fi
		if [ ! -d "$Samplefolder/rRNA" ] ; then
                        mkdir "$Samplefolder/rRNA"
                fi
	fi
	if [ ! -d "$Countsfolder" ] ; then
		mkdir -p "$Countsfolder/RAW"
		mkdir "$Countsfolder/Simple"
		mkdir "$Countsfolder/tRNA_RAW"
		mkdir "$Countsfolder/tRNA_anti"
	else
		if [ ! -d "$Countsfolder/RAW" ] ; then
			mkdir "$Countsfolder/RAW"
		fi
		if [ ! -d "$Countsfolder/Simple" ] ; then
			mkdir "$Countsfolder/Simple"
		fi
		if [ ! -d "$Countsfolder/tRNA_RAW" ] ; then
			mkdir "$Countsfolder/tRNA_RAW"
		fi
		if [ ! -d "$Countsfolder/tRNA_anti" ] ; then
			mkdir "$Countsfolder/tRNA_anti"
		fi
	fi
	if [ ! -d "$Trimfolder" ] ; then
		mkdir "$Trimfolder"
	fi
fi

#2 Start adapt trimming,
file="$RAWfolder/$Sample*.fastq.gz"
file1=$(echo $file|cut -f 1 -d " ")
file2=$(echo $file|cut -f 2 -d " ")
num_of_file=$(ls $file|wc -l)
if [ $num_of_file != 2 ];then echo "More than 2 files share the same sample name, make sure the name is uniq and you have PE reads"; exit 1;fi
trimed1="$Trimfolder/$Sample.1.fq.gz"
trimed2="$Trimfolder/$Sample.2.fq.gz"
#cutadapt -m 15 -O 5 -n 3 -q 20 -b AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -B GATCGTCGGACTGTAGAACTCTGAACGTGTAGA -o $trimed1 -p $trimed2 $file1 $file2

#3 Hisat2 mapping, Pass1
hisat2 -p 24 -k 10 --no-mixed --no-discordant --known-splicesite-infile $Splicesite --novel-splicesite-outfile "$Samplefolder/Hisat/novelsite.txt" -x $Ref -1 $trimed1 -2 $trimed2 |samtools view -bS - > "$Samplefolder/Hisat/hisat.bam"

#4 Process Tophat output
cd "$Samplefolder/Hisat/"
samtools view -H hisat.bam > header.sam
# remove non-human and non-concordant reads
samtools view -F4 hisat.bam |grep -w "NH:i:1" | awk '{CIGAR1=$6;L1=$0;chr1=$3;getline;CIGAR2=$6;L2=$0;chr2=$3;if \
((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)) print L1"\n"L2}' > uniq.sam
samtools view -F4 hisat.bam |grep -v -w "NH:i:1" | awk '{CIGAR1=$6;L1=$0;chr1=$3;getline;CIGAR2=$6;L2=$0;chr2=$3;if \
((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)) print L1"\n"L2}' > multi.sam

#5 Preparing reads for Bowtie2 mapping
cd "$Samplefolder/Bowtie/"
samtools view -@ 24 -bf4 "$Samplefolder/Hisat/hisat.bam" > unmapped.bam
samtools bam2fq -1 unmapped.1.fq -2 unmapped.2.fq -s unmapped.unpaired.fq unmapped.bam
rm unmapped.bam
gzip -f unmapped.1.fq
gzip -f unmapped.2.fq
echo 'Finished HISAT2: ' $Sample

#6 Bowtie2 mapping (Pass2)
bowtie2 --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 -p 24 -k 10 \
	--no-mixed --no-discordant -x $Bowtie -1 unmapped.1.fq.gz -2 unmapped.2.fq.gz \
| samtools view -bS - > bowtie2.bam

#7 Process Bowtie2 output
samtools view -bf4 -@ 24 bowtie2.bam > unmapped.bam
# remove non-human and non-concordant reads, and those with more than 10 bases being soft-clipped in either end from either read
samtools view -F4 -q255 bowtie2.bam |awk '{CIGAR1=$6;L1=$0;chr1=$3;getline;CIGAR2=$6;L2=$0;chr2=$3;if \
((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)) print L1"\n"L2}' > uniq.sam
samtools view -F4 bowtie2.bam |awk '{if ($5<255) print }' | awk '{CIGAR1=$6;L1=$0;chr1=$3;getline;CIGAR2=$6;L2=$0;chr2=$3;if \
((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)) print L1"\n"L2}' > multi.sam
echo 'Finished BOWTIE2: ' $Sample

#8 Combine Pass1 and Pass2, process the protein reads, sense protein reads, and calculate RNAseqmatrix, intersect tRNA reads
cd "$Samplefolder/Combined/"
cat ../Hisat/header.sam ../Hisat/multi.sam ../Bowtie/multi.sam | samtools view -b > multi.bam
python ~/tgirt_benchmark/pipelines/genome_mapping/yaojun_pipeline/reduce_multi_reads.py  \
	--infile multi.bam \
	--outfile multi_filtered.bam \
	--bam_in --bam_out 
samtools view multi_filtered.bam \
	| cat ../Hisat/header.sam ../Hisat/uniq.sam ../Bowtie/uniq.sam - \
	| samtools view -@24 -bS -	\
	| samtools sort -n -@ 24 -O bam -T temp \
	> primary.bam
echo "Finished correcting multimple mapped reads: " $Sample

bedtools pairtobed -s -f 0.01 -abam primary.bam -b $REF/GRCh38/Bed_for_counts_only/tRNA.bed > ../tRNA/tRNA_primary.bam
bedtools pairtobed -s -f 0.01 -abam primary.bam -b $REF/GRCh38/Bed_for_counts_only/rRNA_for_bam_filter.bed > ../rRNA/rRNA_primary.bam
bedtools pairtobed -s -f 0.01 -abam primary.bam -b $REF/GRCh38/Bed_for_counts_only/sncRNA_no_tRNA.bed > sncRNA.bam
bedtools pairtobed -s -f 0.01 -type neither -abam primary.bam -b $REF/GRCh38/Bed_for_counts_only/sncRNA_rRNA_for_bam_filter.bed > primary_no_sncRNA_tRNA_rRNA.bam
bedtools pairtobed -abam primary_no_sncRNA_tRNA_rRNA.bam -b $REF/GRCh38/Bed_for_counts_only/protein.bed > protein.bam
mkdir -p temp
cd temp
samtools view -@ 24 -bF64 ../protein.bam > R2.bam
samtools view -@ 24 -bf64 ../protein.bam > R1.bam
bedtools intersect -s -wa -a R1.bam -b $REF/GRCh38/Bed_for_counts_only/protein.bed > R1_1.bam
bedtools intersect -S -wa -a R2.bam -b $REF/GRCh38/Bed_for_counts_only/protein.bed > R2_1.bam
bedtools intersect -s -bed -v -f 0.01 -wa -a R1_1.bam -b $REF/GRCh38/Bed_for_counts_only/sncRNA_x_protein.bed |cut -f 4|cut -f 1 -d "/" > R1.id
bedtools intersect -S -bed -v -f 0.01 -wa -a R2_1.bam -b $REF/GRCh38/Bed_for_counts_only/sncRNA_x_protein.bed |cut -f 4|cut -f 1 -d "/" > R2.id
cat R1.id R2.id |sort -u > id.txt
#picard FilterSamReads INPUT=../protein.bam FILTER=includeReadList READ_LIST_FILE=id.txt OUTPUT=../protein.sense.bam WRITE_READS_FILES=false SORT_ORDER=unsorted
cd ..
rm -r temp
#samtools sort -@ 24 -O bam -T temp protein.sense.bam > temp.bam
#mv temp.bam protein.sense.bam
echo 'Made protein bam and split type: ' $Sample

bedtools bamtobed -mate1 -bedpe -i sncRNA.bam > sncRNA.bedpe
awk '{FS="\t"; OFS="\t"; if ($2>$5) $2=$5; if ($3<$6) $3=$6; print $1,$2,$3,$7,0,$9}' sncRNA.bedpe > sncRNA.bed
awk '{print $3-$2}' sncRNA.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > sncRNA_span.txt
bedtools bamtobed -mate1 -bedpe -i primary_no_sncRNA_tRNA_rRNA.bam > primary_no_sRNAs.bedpe
awk '{FS="\t"; OFS="\t"; if ($2>$5) $2=$5; if ($3<$6) $3=$6; print $1,$2,$3,$7,0,$9}' primary_no_sRNAs.bedpe > primary_no_sRNAs.bed
awk '{print $3-$2}' primary_no_sRNAs.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > primary_no_sRNAs_span.txt
echo 'Converted to bed' $Sample

#9 tRNA reads process
cd "$Samplefolder/tRNA/"
samtools bam2fq -1 tRNA.1.fq -2 tRNA.2.fq -s tRNA.unpaired.fq tRNA_primary.bam
gzip -f tRNA.1.fq
gzip -f tRNA.2.fq
bowtie2 -p 24 --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 --norc --no-mixed --no-discordant -x $REF/GRCh38/tRNA/tRNA -1 tRNA.1.fq.gz -2 tRNA.2.fq.gz | samtools view -bS - > tRNA_remap.bam
samtools view -bF4 tRNA_remap.bam | bedtools bamtobed -mate1 -bedpe -i - > tRNA.bedpe
awk '{FS="\t"; OFS="\t"; if ($2>$5) $2=$5; if ($3<$6) $3=$6; print $1,$2,$3,$7,0,$9}' tRNA.bedpe > tRNA.bed
awk '{print $3-$2}' tRNA.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > tRNA_span.txt
awk '{print $2+1}' tRNA.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > tRNA_start.txt
awk '{print $3}' tRNA.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > tRNA_end.txt
echo "Mapped tRNA reads: " $Sample

#10 Generate tRNA counts
samtools view tRNA_remap.bam|cut -f 3|grep -v "*"|sort|uniq -c|awk '{print $2"\t"$1/2}' > "$Countsfolder/tRNA_RAW/$Sample.tRNA"
tRNA_file="$Countsfolder/tRNA_RAW/$Sample.tRNA"
tRNA_anti="$Countsfolder/tRNA_anti/$Sample.anti"
rm -f $tRNA_anti
anticodon="AlaAGC AlaCGC AlaTGC GlyGCC GlyCCC GlyTCC ProAGG ProCGG ProTGG ThrAGT ThrCGT ThrTGT\
        ValAAC ValCAC ValTAC PheGAA AsnATT AsnGTT LysCTT LysTTT AspGTC GluCTC GluTTC HisGTG\
        GlnCTG GlnTTG SerAGA SerCGA SerTGA SerACT SerGCT ArgACG ArgCCG ArgTCG ArgCCT ArgTCT LeuAAG\
        LeuCAG LeuTAG LeuCAA LeuTAA IleAAT IleGAT IleTAT TyrATA TyrGTA CysGCA CysACA TrpCCA\
        Undet??? SupCTA SupTTA"
for i in `echo $anticodon`
        do
        sum=0
        grep $i $tRNA_file | awk '{sum+=$2}END{if (sum>0) print "'"$i"'\t"sum; else print "'"$i"'\t"0}' 
done > $tRNA_anti
sum=0
grep 'SeCTCA\|SeC(e)TCA' $tRNA_file|awk '{sum+=$2}END{if (sum>0) print "SelCysTCA\t"sum; else print "SelCysTCA\t"0}' >> $tRNA_anti
sum=0
grep '10-MetCAT\|16.trna20-MetCAT\|75-MetCAT\|97-MetCAT\|164-MetCAT\|162-MetCAT\|27-MetCAT\|21-MetCAT\|22-MetCAT\|92-MetCAT' $tRNA_file|awk '{sum+=$2}END{if (sum>0) print "MetCAT\t"sum; else print "MetCAT\t"0}' >> $tRNA_anti
sum=0
grep '32-MetCAT\|17.trna20-MetCAT\|2-MetCAT\|171-MetCAT\|169-MetCAT\|150-MetCAT\|142-MetCAT\|129-MetCAT\|61-MetCAT\|1-MetCAT' $tRNA_file|awk '{sum+=$2}END{if (sum>0) print "iMetCAT\t"sum; else print "iMetCAT\t"0}' >> $tRNA_anti
echo 'Finished counting tRNA:' $Sample

#12 rRNA reads process
cd "$Samplefolder/rRNA/"
samtools bam2fq -1 rRNA.1.fq -2 rRNA.2.fq -s rRNA.unpaired.fq rRNA_primary.bam
gzip -f rRNA.1.fq
gzip -f rRNA.2.fq
bowtie2 -p 24 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 --no-mixed --no-discordant -x $REF/GRCh38/Plasma_ref/rDNA -1 rRNA.1.fq.gz -2 rRNA.2.fq.gz | samtools view -bS - > rRNA_remap.bam
samtools view -bF4 rRNA_remap.bam | bedtools bamtobed -mate1 -bedpe -i - > rRNA.bedpe
awk '{FS="\t"; OFS="\t"; if ($2>$5) $2=$5; if ($3<$6) $3=$6; print $1,$2,$3,$7,0,$9}' rRNA.bedpe > rRNA.bed
awk '{print $3-$2}' rRNA.bed |sort|uniq -c| sort -k 2n|awk '{print $2,"\t",$1}' > rRNA_span.txt
bedtools coverage -s -counts -F 0.1 -a $REF/GRCh38/Bed_for_counts_only/rDNA.bed -b rRNA.bed > rRNA.counts
echo 'Finished counting rRNA:' $Sample

#13 Generate gene counts file
cd "$Samplefolder/Combined/"
bedtools coverage -s -counts -F 0.1 -a $REF/GRCh38/Bed_for_counts_only/sncRNA_no_tRNA.bed -b sncRNA.bed > sncRNA.counts
bedtools coverage -s -counts -F 0.1 -a $REF/GRCh38/Bed_for_counts_only/genes_no_sncRNA_rRNA_tRNA.bed -b primary_no_sRNAs.bed  > non_sRNAs.counts
bedtools coverage -s -counts -F 0.1 -a $REF/RNASeqConsortium/ercc/ercc.bed -b primary_no_sRNAs.bed > ercc.counts

cd "$Countsfolder"
cat ../$Sample/Combined/non_sRNAs.counts ../$Sample/Combined/sncRNA.counts ../$Sample/rRNA/rRNA.counts ercc.counts > RAW/$Sample.counts
awk '{print $8,$1,$2,$3,$6,$4,$7,$9}' OFS=\\t FS=\\t RAW/$Sample.counts > Simple/$Sample.s.counts
cd Simple
sed -i 's/3prime_overlapping_ncrna/other_ncRNA/g' $Sample.s.counts
sed -i 's/processed_transcript/other_ncRNA/g' $Sample.s.counts
sed -i 's/sense_intronic/other_ncRNA/g' $Sample.s.counts
sed -i 's/sense_overlapping/other_ncRNA/g' $Sample.s.counts
sed -i 's/IG_.*gene/IG/g' $Sample.s.counts
sed -i 's/TR_.*gene/TR/g' $Sample.s.counts
sed -i 's/Mt_.RNA/Mt/g' $Sample.s.counts
sed -i 's/Mt_protein_coding/Mt/g' $Sample.s.counts
sed -i 's/trans.*pseudogene/pseudogene/g' $Sample.s.counts
sed -i 's/un.*pseudogene/pseudogene/g' $Sample.s.counts
sed -i 's/p.*pseudogene/pseudogene/g' $Sample.s.counts
awk '{if (($7=="misc_RNA")&&($6~/Y_RNA|RNY/)) print $1,$2,$3,$4,$5,$6,"YRNA",$8;else print $1,$2,$3,$4,$5,$6,$7,$8}' FS=\\t OFS=\\t $Sample.s.counts > temp.txt
mv temp.txt $Sample.s.counts
awk '{if (($7=="misc_RNA")&&($6~/VTRNA|Vault/)) print $1,$2,$3,$4,$5,$6,"VTRNA",$8;else print $1,$2,$3,$4,$5,$6,$7,$8}' FS=\\t OFS=\\t $Sample.s.counts > temp.txt
mv temp.txt $Sample.s.counts
awk '{if (($7=="misc_RNA")&&($6~/7SL|SRP/)) print $1,$2,$3,$4,$5,$6,"7SL",$8;else print $1,$2,$3,$4,$5,$6,$7,$8}' FS=\\t OFS=\\t $Sample.s.counts > temp.txt
mv temp.txt $Sample.s.counts
awk '{if (($7=="misc_RNA")&&($6~/7SK/)) print $1,$2,$3,$4,$5,$6,"7SK",$8;else print $1,$2,$3,$4,$5,$6,$7,$8}' FS=\\t OFS=\\t $Sample.s.counts > temp.txt
mv temp.txt $Sample.s.counts
sort -k 7,7d -k 1,1d $Sample.s.counts > $Sample.ver2.s.counts
mv $Sample.ver2.s.counts $Sample.s.counts
echo 'Finished counting all gene:' $Sample


cd "$Samplefolder/Hisat/"
rm *.sam
cd "$Samplefolder/Bowtie/"
rm *.sam
cd "$Samplefolder/Combined/"
#rm *.sam
#picard CollectRnaSeqMetrics REF_FLAT=$REF/GRCh38/hg38_rDNA/hg38.refflat STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=20 INPUT=protein.sense.bam OUTPUT=RnaSeq.Metrics REFERENCE_SEQUENCE="$Bowtie.fa" ASSUME_SORTED=false
