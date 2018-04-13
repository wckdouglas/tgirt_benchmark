#!/usr/bin/env Rscript

#library("BiocParallel")
#register(MulticoreParam(12))
library(DESeq2)
library(cowplot)
library(feather)
library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(purrr)

gene_df <- read_tsv('/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv') %>%
    dplyr::rename(id = gene_id) %>%
    dplyr::select(name, id, type) %>% 
    mutate(id = str_replace(id, '_gene$','')) %>%
    distinct()

#subset data table to retain needed columns
selectSample <- function(d, pattern) {
    d <- d %>% dplyr::select(grep(str_c('id|',pattern),names(.)))
    d <- data.frame(d)
    row.names(d) <- d[,1]
    return (d[,-1])
}


# fitting deseq to table using input mix comparison
fitDESeq <-  function(comparison,df){
    df <- selectSample(df, comparison)
    col_data <- data.frame(samplenames = names(df)) %>% 
        mutate(mix = str_sub(samplenames,8,8)) %>% 
        mutate(mix =  factor(mix,levels = rev(unique(mix)))) 
    row.names(col_data) <- col_data$samplenames
    de <- DESeqDataSetFromMatrix(countData = df,
                           colData = col_data,
                           design = ~mix) %>%
        DESeq() %>%
        results() %>%
        data.frame  %>%
        rownames_to_column(var = "id") %>%
        mutate(comparison = str_replace(comparison,'\\|',' vs ')) %>%
        tbl_df
    return(de)
}

read_table_func <- function(tablename){
    tablename %>%
        read_tsv()  %>%
        mutate(id = str_replace(id, '_gene$','')) %>%
        mutate(id = str_replace(id, '[0-9]+-[0-9]+$','')) %>%
        mutate(id = str_replace(id,'\\([+-]\\)','')) %>%
        left_join(gene_df) %>% 
        mutate(id = ifelse(type %in% c('Mt','rRNA'),name, id)) %>%
        dplyr::select(-name, -type) %>%
        group_by(id) %>% 
        summarize_all(sum) %>%
        ungroup() %>%
        tbl_df %>%
        return()
}

project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'
out_path <- file.path(project_path, '/DEgenes')
read_table_and_DESeq <- function(tablename, map_type){
    out_table <- str_c(out_path,'/', map_type,'.feather')
    tablename %>%
    read_table_func() %>%
    map_df(c('A|B','C|D'), fitDESeq, .) %>%
    mutate(map_type = map_type)  %>%
    write_feather(out_table)
    message('Written: ', out_table)
    return(0)
}
    
#table_names <- c('/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/Counts/RAW/combined_gene_count.tsv',
table_names <- c(file.path(project_path,'genome_mapping/tgirt_map/Counts/RAW/combined_gene_count.tsv'),
                 file.path(project_path,'genome_mapping/conventional/counts/feature_counts.tsv'))
map_types <- c('Customized_pipeline','Conventional_pipeline')
result <- map2(table_names, map_types, read_table_and_DESeq)


