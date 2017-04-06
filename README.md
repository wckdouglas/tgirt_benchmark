# TGIRT-seq benchmarking #

This is the repository for storing scripts attempting to compare pipelines for TGIRT-seq gene quantifications.

We will use MAQC samples for test data.

---

## Download data ##

```
python download_data/download.py
```

and run all the command it generates

## Make reference ##

### make transcripts reference ###

1. download Ensembl transcripts
2. download ERCC transcripts
3. get rRNA fasta file from NCBI
4. merge with customized tRNA fasta
5. Make table matching transcripts to genes
6. bowtie2 index [rRNA and tRNA]

```
cd download_data
bash make_transcript_ref.sh
```

### make genome reference ###

1. downlaod ensembl genes from biomart
2. Make bed file with rRNA, tRNA and ERCC  and all ensembl genes
3. Make SAF file for **featureCounts**
4. split bed file for counts
5. download Ensembl hg38 genome
6. spike in rRNA and ERCC fasta
7. make **HISAT2** index from the reference

```
bash make_genome_ref.sh
```

---

## Aligning and counts ##

### Run Kallisto ###


General command is as followed:

```
kallisto quant \
	-i ${INDEX} -o ${COUNT_PATH}/${SAMPLENAME} \
	--fr-stranded  --threads=${THREADS}\
	${R1} ${R2}
```

```
cd pipelines/alignment_free/kallisto
bash build_index.sh
bash kallisto.sh > command.sh
bash command.sh
```

### Run Salmon ###

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
cd pipelines/alignment_free/salmon
bash build_index.sh
bash salmon.sh > command.sh
bas command.sh
```

### Run genome mapping with customized maps

```
bash prepair_jun_Commands.sh > command.sh
bash command.sh
```
