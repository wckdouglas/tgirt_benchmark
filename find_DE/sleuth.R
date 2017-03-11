#!/usr/bin/env Rscript


options(mc.cores = 12)
library(sleuth)
library(readr)
library(stringr)
library(purrr)
library(dplyr)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/GRCh38/transcripts.tsv' %>%
    read_tsv()  %>%
    dplyr::rename(target_id=t_id) %>%
    tbl_df

gene_name <- gene_file %>%
    select(-target_id) %>%
    unique()

# make sample dataframe
project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking/alignment_free/countFiles'
count_files_df <- list.files(project_path) %>%
    data.frame(sample_id = .) %>%
    filter(!grepl('.tsv$', sample_id)) %>%
    mutate(sample_mix = str_sub(sample_id, 8,8)) %>%
    mutate(sample_number = str_sub(sample_id, 10,10)) %>%
    mutate(path = str_c(project_path, sample_id, sep='/')) %>%
    dplyr::rename(sample=sample_id) %>%
    tbl_df

run_sleuth_DE <- function(sample_regex){
    mixes <- str_split(sample_regex,'\\|')[[1]]
    so <- count_files_df %>%
        filter(grepl(sample_regex, sample_mix)) %>%
        sleuth_prep( ~ sample_mix, 
                     target_mapping = gene_file, 
                     aggregation_column = 'gene_id') %>%
        sleuth_fit() %>%
        sleuth_fit(~1, 'reduced') %>%
        sleuth_lrt('reduced', 'full') %>%
        sleuth_results('reduced:full', test_type = 'lrt')
        dplyr::rename(gene_id=target_id) %>%
        mutate(sample_base = mixes[1]) %>% 
        mutate(sample_test = mixes[2]) %>%
        inner_join(gene_name , by = 'gene_id')
    return(so)
}

sleuth_AB <- run_sleuth_DE('A|B')
sleuth_CD <- run_sleuth_DE('C|D')
