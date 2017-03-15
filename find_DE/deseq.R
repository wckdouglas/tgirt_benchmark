#!/usr/bin/env Rscript

library(DESeq2)
library(cowplot)
library(feather)
library(readr)
library(stringr)
library(tibble)
library(dplyr)

selectSample <- function(d, pattern) {
    d <- d %>% select(grep(str_c('id|',pattern),names(.)))
    d <- data.frame(d)
    row.names(d) <- d[,1]
    return (d[,-1])
}

fitDESeq <-  function(df){
    sdf <- data_frame(samplenames = names(df)) %>%
        separate(samplenames,c('sample','biosample','replicate'),sep='_') %>%
        data.frame()
    col_data <- data.frame(samplenames = names(df)) %>% 
        mutate(annotation = str_sub(samplenames,8,8)) %>% 
        mutate(annotation =  factor(annotation,levels = rev(unique(annotation)))) 
    row.names(col_data) <- col_data$samplenames
    de <- DESeqDataSetFromMatrix(countData = df,
                           colData = col_data,
                           design = ~annotation) %>%
        DESeq() %>%
        results() %>%
        data.frame  %>%
        select(log2FoldChange,padj,baseMean) %>%
        setNames(paste(c('logFC','adj.P.Val','AveExpr'), 
                 paste(unique(sdf$biosample),collapse = ''),sep='_')) %>%
        rownames_to_column(var = "id") %>%
        tbl_df
    return(de)
}

rename_ryan <- function(column_name){
    if (grepl('ref',column_name)){
        sample_mix <- str_sub(column_name,4,4)
        sample_id <- str_sub(column_name,5,5)
        return(str_c('Sample',sample_mix,sample_id,sep='_'))
    }else{
        return(column_name)
    }
}

# read files
project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking'
genome_df <- project_path %>%
    str_c('/genome_mapping/RAW/combined_gene_count.tsv') %>%
    read_tsv()

ryan_df <- '/stor/work/Lambowitz/Data/archived_work/TGIRT_ERCC_project/result/countTables' %>%
    str_c('countsData.tsv', sep='/') %>%
    read_tsv() %>%
    select(grep('ref|id|type|name', names(.)))  %>%
    set_names(sapply(names(.),rename_ryan)) %>%
    mutate(id = ifelse(type=='tRNA', str_replace(id,'[0-9]+',''),id)) %>%
    mutate(name = ifelse(type=='tRNA', str_replace(name,'[0-9]+',''),name)) %>%
    group_by(id,type,name) %>%
    summarise_all(sum) %>%
    ungroup() %>%
    tbl_df
    
# ================ Make sample table ===================
g_AB <- selectSample(genome_df, 'A|B')
g_CD <- selectSample(genome_df, 'C|D')
ryan_AB <- selectSample(ryan_df,'A|B')
ryan_CD <- selectSample(ryan_df,'C|D')

deseq_fc_df <- inner_join(fitDESeq(g_AB), fitDESeq(g_CD)) %>%
    mutate(map_type = 'W/ multimap') 
ryan_deseq_fc_df <- inner_join(fitDESeq(ryan_AB), fitDESeq(ryan_CD)) %>%
    mutate(map_type = 'W/o multimap') 
df <- deseq_fc_df %>%
    rbind(ryan_deseq_fc_df)
out_table <- str_c(project_path,'/genome_mapping/deseq_genome.feather',sep='/')
write_feather(df, out_table)