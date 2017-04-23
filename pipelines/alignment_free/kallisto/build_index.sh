cd $REF/benchmarking/human_transcriptome
KALLISTO_REF=transcriptome_kallisto
mkdir $KALLISTO_REF
cd $KALLISTO_REF
kallisto index -i transcriptome_kallisto ../whole_transcriptome.fa
