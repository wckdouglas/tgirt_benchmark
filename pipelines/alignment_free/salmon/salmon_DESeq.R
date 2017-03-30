#!/usr/bin/env Rscript

library(readr)
library(stringr)
library(purrr)
library(dplyr)
library(tximport)
library(feather)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/GRCh38/transcripts.tsv' %>%
    read_tsv()  %>%
    dplyr::rename(target_id=t_id) %>%
    tbl_df

tx2gene <- gene_file %>%
    select(target_id, gene_id) %>%
    set_names(c('TXNAME','GENEID'))


# tximport salmon abundance to gene count
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/alignment_free/salmon'
salmon_files <-  list.files(project_path, pattern = '[1-3]$',full.names = T) %>%
    str_c('/quant.sf')
salmon_df <- tximport(kallisto_files, 
                        type = "salmon", 
                        tx2gene = tx2gene, 
                        reader = read_tsv)