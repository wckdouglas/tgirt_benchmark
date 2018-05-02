cd $REF/benchmarking/human_transcriptome
for kmer in 11 15 21
do
	salmon index -i transcript_salmon_$kmer \
		--perfectHash \
		--threads 24 \
		--transcripts whole_transcriptome.fa \
		--kmerLen $kmer
done

salmon index -i transcript_salmon \
    --perfectHash \
    --threads 24 \
    --transcripts whole_transcriptome.fa \
    --kmerLen 31
