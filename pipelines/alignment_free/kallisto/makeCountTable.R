#!/usr/bin/env Rscript

library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

readTPM <- function(samplename, kallisto_dir) {
    filename <- str_c(kallisto_dir,'/',samplename,'/abundance.tsv')
    df <- read_tsv(filename) %>%
        select(tpm, target_id) %>%
        mutate(sample_id = str_replace_all(samplename,'-','_')) 
    message('Read: ',samplename)
    return(df)
}


kallisto_dir <- '/stor/scratch/Lambowitz/cdw2854/bench_marking/alignment_free/countFiles'
count_table <-  str_c(kallisto_dir,'/alignment_free.tsv')
df <- list.files(kallisto_dir) %>%
	.[!grepl('.tsv$|.feather$',.)] %>%
    map(readTPM, kallisto_dir) %>%
    reduce(rbind) %>%
    spread(sample_id, tpm) %>%
    write_tsv(count_table)
message('Written ', count_table)
