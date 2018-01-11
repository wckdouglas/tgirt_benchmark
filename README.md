# TGIRT-seq benchmarking #

This is the repository for storing scripts attempting to compare pipelines for TGIRT-seq gene quantifications.

We used [MAQC samples](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=2017824) as test data.

---

## Download data ##

```
python download_data/download.py
```

and run all the command it generates

## Make reference ##

1. download Ensembl transcripts
2. download ERCC transcripts
3. get rRNA fasta file from NCBI
4. merge with customized tRNA fasta
5. Make table matching transcripts to genes
6. bowtie2 index [rRNA and tRNA]
7. downlaod ensembl genes from biomart
8. Make bed file with rRNA, tRNA and ERCC  and all ensembl genes
9. Make SAF file for **featureCounts**
10. split bed file for counts
11. download Ensembl hg38 genome
12. spike in rRNA and ERCC fasta
13. make **HISAT2** index from the reference


```
cd download_data
bash make_reference.sh
```

---

## Aligning and/or quantification ##

All path should be changed accordingly if running on other machines.

### Alignment-free pipelines ###

#### Run Kallisto ####

All scripts are under: ```pipelines/alignment_free/kallisto/```

General command is as followed:

```
kallisto quant \
	-i ${INDEX} -o ${COUNT_PATH}/${SAMPLENAME} \
	--fr-stranded  --threads=${THREADS}\
	${R1} ${R2}
```

```
bash build_index.sh
bash kallisto.sh > command.sh
bash command.sh
Rscript kallisto_DESeq.R          #Tximport then DESeq2
```

#### Run Salmon ####

All scripts are under: ```pipelines/alignment_free/salmon/```

General command is as followed:

```
salmon quant \
        --seqBias --gcBias \
        --index $INDEX --libType ISF \
        --writeMappings  \
        --threads=${THREADS} \
        --auxDir $RESULTPATH/${SAMPLENAME} \
        --numBootstraps 100 \
    	--mates1 ${R1} --mates2 ${R2} \
		--output $RESULTPATH/${SAMPLENAME} \
    \| samtools view -b@ $THREADS \
    \> ${RESULTPATH}/${SAMPLENAME}.bam
```


```
bash build_index.sh
bash salmon.sh > command.sh
bas command.sh
Rscript salmon_DESeq.R    #Tximport then DESeq2
```

### Alignment-based pipelines ###

#### Run hisat2+featureCounts ####

All scripts are under: ```pipelines/genome_mapping/conventional/```

```
bash hisat2.sh > command.sh
bash command.sh
python feature_counts.py 
python clean_output.py    # for cleaning the featureCounts output for DESeq2
```


#### Run genome mapping with TGIRT-map ####

All scripts are under: ```pipelines/genome_mapping/tgirt_map_pipeline/```

The pipeline can be downloaded [here](https://github.com/wckdouglas/tgirt_map)

```
prepair_command.sh > command.sh
bash command.sh
python makeCountTable.py       #for making a gene count table
```

#### Fold-change analysis ####

```
Rscript pipelines/genome_mapping/deseq.R
```


## Plottings ##

Analysis and plottings are mostly done in R scripts under ```find_DE/```.
