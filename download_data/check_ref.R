#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)

bed_path <- '/stor/work/Lambowitz/ref/RNASeqConsortium'
gene_bed <- bed_path %>%
    str_c('/genes.bed') %>%
    read_tsv(col_names=c('chrom','start','end','gene_name','score','strand','bio_type','gene_id')) %>%
    tbl_df

df <- read_tsv('/stor/work/Lambowitz/ref/human_transcriptome/transcripts.tsv') %>%
     right_join(gene_bed)

df %>% filter(is.na(name))  %>% select(-chrom,-start,-end, -score, -strand) %>% View
