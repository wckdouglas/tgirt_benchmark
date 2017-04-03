#!/usr/bin/env Rscript

library("BiocParallel")
register(MulticoreParam(12))
library(readr)
library(stringr)
library(purrr)
library(dplyr)
library(tximport)
library(feather)
library(DESeq2)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/GRCh38/transcripts.tsv' %>%
    read_tsv()  %>%
    dplyr::rename(target_id=t_id) %>%
    tbl_df

tx2gene <- gene_file %>%
    select(target_id, gene_id) %>%
    set_names(c('TXNAME','GENEID'))

# make sample file and annotations
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/alignment_free/kallisto/countFiles'
kallisto_files_df <-  list.files(project_path, pattern = '[1-3]$') %>%
    data.frame(samplename=.) %>%
    mutate(filename = str_c(project_path,samplename,'abundance.tsv',sep='/'))%>%
    mutate(samplename = str_replace_all(samplename,'-','_')) %>%
    mutate(mix = str_sub(samplename, 8,8)) %>%
    mutate(sample_id = str_sub(samplename, 11, 10)) %>%
    tbl_df
    

# fit deseq to selected samples (A vs B ['AB']and C vs D ['CD'])
fit_DESeq <- function(sample_comparison){
    kallisto_files_df <- filter(kallisto_files_df, grepl(sample_comparison,mix))
    kallisto_files <- kallisto_files_df$filename
    names(kallisto_files) <- kallisto_files_df$samplename

    # condition data frame for deseq2
    cond_df <- kallisto_files_df %>%
        select(mix, samplename) %>%
        mutate(mix =  factor(mix,levels = rev(unique(mix))))  %>%
        data.frame()
    rownames(cond_df) = cond_df$samplename

    # tximport kallisto abundance to gene count
    kallisto_df <- tximport(kallisto_files, 
                        type = "kallisto", 
                        tx2gene = tx2gene, 
                        reader = read_tsv)

    # run deseq2 on tximport table
    dds <- DESeqDataSetFromTximport(kallisto_df, cond_df, ~mix) %>%
        DESeq() %>%
        results() %>%
        data.frame %>%
        rownames_to_column(var = "id") %>%
        mutate(comparison = str_replace(sample_comparison,'\\|',' vs ')) %>%
        tbl_df
    return(dds)
}

out_file_name = file.path(dirname(project_path),'kallisto_DESeq.feather')
kallisto_df <- map(c('A|B','C|D'), fit_DESeq)  %>%
    purrr:::reduce(rbind) %>%
    write_feather(out_file_name)
message('Written: ', out_file_name)
