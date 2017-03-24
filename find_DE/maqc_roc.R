#!/usr/bin/evn Rscript

library(feather)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

taqman <- '/stor/scratch/Lambowitz/cdw2854/bench_marking/maqc/taqman_fc_table.feather' %>%
    read_feather() %>%
    rename(id = ensembl_gene_id) %>%
    mutate(real_FC = ifelse(abs(logFC_AB)>0.5, 'yes','no')) %>%
    select(-logFC_AB)

project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking'
alignment_free_df <- project_path %>%
    str_c('/alignment_free/countFiles/sleuth_results.feather') %>%
    read_feather() %>%
    mutate(logFC = -b) %>%
    dplyr::rename(id = gene_id) %>%
    select(id,sample_base,sample_test, name,type,qval, mean_obs, logFC, pval) %>%
    gather(variable, value, -id,-name,-type,-sample_base,-sample_test) %>%
    mutate(variable = case_when(
        .$variable == 'qval' ~ 'adj.P.Val',
        .$variable == 'mean_obs' ~ 'AveExpr',
        .$variable == 'logFC' ~ 'logFC',
        .$variable == 'pval' ~ 'P.Val'
    ))  %>%
    mutate(variable = str_c(variable,'_', sample_base, sample_test)) %>%
    select(-sample_base,-sample_test) %>%
    spread(variable, value)  %>%
    select(-name,-type) %>%
    mutate(map_type = 'Kallisto') %>%
    tbl_df

genome_df <- project_path %>%
    str_c('/genome_mapping/pipeline7_counts/deseq_genome.feather') %>%
    read_feather() %>%
    tbl_df
df <- rbind(genome_df, alignment_free_df) %>%
    inner_join(taqman)