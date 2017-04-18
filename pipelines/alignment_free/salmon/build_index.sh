cd $REF/human_transcriptome
salmon index -i transcript_salmon --perfectHash --threads 24 --transcripts whole_transcriptome.fa
