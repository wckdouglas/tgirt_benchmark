#!/bin/bash

DATAPATH=${SCRATCH}/bench_marking/genome_mapping/tgirt_map/Trim
REF_PATH=${SCRATCH}/ref/benchmarking_new
INDEX_PATH=${REF_PATH}/benchmarking/human_transcriptome
TRANSCRIPTOME_INDEX=$INDEX_PATH/whole_transcriptome
TRANSCRIPTOME_FASTA=${TRANSCRIPTOME_INDEX}.fa
RESULTPATH=${SCRATCH}/bench_marking/alignment_free/salmon_aligned
BAM_PATH=$RESULTPATH/bam_files
THREADS=24
mkdir -p $RESULTPATH $BAM_PATH

for R1 in ${DATAPATH}/*.1.fq.gz
do
    R2=${R1/.1./.2.}
    SAMPLENAME=$(basename ${R1%.1.fq.gz})
    BAM_FILE=$BAM_PATH/${SAMPLENAME}.bam
    echo bowtie2 --threads $THREADS \
            -x $TRANSCRIPTOME_INDEX \
            --score-min G,1,10 \
            -k 30 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 \
            -1 ${R1} -2 ${R2} \
            --no-mixed --no-discordant --dovetail \
            --fr --local  \
        \| samtools view -b@ $THREADS \
        \> $BAM_FILE \
        \; time salmon quant \
            --seqBias --gcBias \
            --libType ISF \
            --targets $TRANSCRIPTOME_FASTA \
            --threads=${THREADS} \
            --auxDir ${RESULTPATH}/${SAMPLENAME} \
            --numBootstraps 100 \
            --useErrorModel  \
            --alignments $BAM_FILE \
            --output ${RESULTPATH}/${SAMPLENAME} 
done

