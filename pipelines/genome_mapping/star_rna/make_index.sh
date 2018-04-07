GENOME_PATH=/stor/work/Lambowitz/ref/benchmarking/GRCH38_genome
STAR_INDEX=$GENOME_PATH/star_index
mkdir -p $STAR_INDEX 

STAR \
    --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir $STAR_INDEX \
    --genomeFastaFiles $GENOME_PATH/reference.fa \
    --sjdbGTFfile $GENOME_PATH/genes.gtf 
