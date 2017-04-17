#!/usr/bin/env Rscript

library(readr)
library(dplyr)
setwd('/stor/work/Lambowitz/ref/human_transcriptome')


gtf_genes <- 'gtf_gene_id' %>%
    read_tsv(col_names = 'gene_id') %>%
    mutate(gtf = 'gtf')

transcript_gene_id <- 'transcripts.tsv' %>%
    read_tsv() 

all_genes <- gtf_genes %>% 
    full_join(transcript_gene_id)


all_genes %>% 
    filter(is.na(name))