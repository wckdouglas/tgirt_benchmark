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

#subset data table to retain needed columns
selectSample <- function(d, pattern) {
    d <- d %>% select(grep(str_c('id|',pattern),names(.)))
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

out_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/DEgenes'
read_table_and_DESeq <- function(tablename, map_type){
    out_table <- str_c(out_path,'/', map_type,'.feather')
    tablename %>%
        read_tsv()  %>%
        tbl_df %>%
        map_df(c('A|B','C|D'), fitDESeq, .) %>%
        mutate(map_type = map_type)  %>%
        write_feather(out_table)
        message('Written: ', out_table)
}
    
table_names <- c('/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/Counts/RAW/combined_gene_count.tsv',
                 '/stor/work/Lambowitz/cdw2854/bench_marking/genome_mapping/Trim/conventional/counts/feature_counts.tsv')
map_types <- c('Customized_pipeline','Conventional_pipeline')
map2(table_names, map_types, read_table_and_DESeq)


