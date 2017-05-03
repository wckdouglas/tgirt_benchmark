#!/usr/bin/env Rscript

#library("BiocParallel")
#register(MulticoreParam(12))
library(readr)
library(stringr)
library(purrr)
library(dplyr)
library(tibble)
library(tximport)
library(feather)
library(DESeq2)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv' %>%
    read_tsv()  %>%
    dplyr::rename(target_id=t_id) %>%
    tbl_df

tx2gene <- gene_file %>%
    select(target_id, gene_id) %>%
    set_names(c('TXNAME','GENEID')) %>%
#    mutate(GENEID = str_replace(GENEID, '[0-9]+\\-[0-9]+','')) %>%
    unique() %>%
    filter(!duplicated(TXNAME))

# make sample file and annotations
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/alignment_free/salmon'
salmon_files_df <-  list.files(project_path, pattern = '[1-3]$') %>%
    data.frame(samplename=.) %>%
    mutate(filename = str_c(project_path,samplename,'quant.sf',sep='/'))%>%
    mutate(samplename = str_replace_all(samplename,'-','_')) %>%
    mutate(mix = str_sub(samplename, 8,8)) %>%
    mutate(sample_id = str_sub(samplename, 11, 10)) %>%
    tbl_df


# fit deseq to selected samples (A vs B ['AB']and C vs D ['CD'])
fit_DESeq <- function(sample_comparison){
    salmon_subset_df <- filter(salmon_files_df, grepl(sample_comparison,mix))
    salmon_files <- salmon_subset_df$filename
    names(salmon_files) <- salmon_subset_df$samplename
    
    # condition data frame for deseq2
    cond_df <- salmon_subset_df %>%
        select(mix, samplename) %>%
        mutate(mix =  factor(mix,levels = rev(unique(mix))))  %>%
        data.frame()
    
    # tximport salmon abundance to gene count
    salmon_df <- tximport(salmon_files, 
                            type = "salmon", 
                            tx2gene = tx2gene, 
                            reader = read_tsv)
    rownames(cond_df) = colnames(salmon_df$counts)
    
    # run deseq2 on tximport table
    dds <- DESeqDataSetFromTximport(salmon_df, cond_df, ~mix) %>%
        DESeq() %>%
        results() %>%
        data.frame %>%
        rownames_to_column(var = "id") %>%
        mutate(comparison = str_replace(sample_comparison,'\\|',' vs ')) %>%
        tbl_df
    return(dds)
}

salmon_df <- map_df(c('A|B','C|D'), fit_DESeq)  %>%
    mutate(map_type = 'Salmon')
out_path <- '/stor/work/Lambowitz/cdw2854/bench_marking/DEgenes'
out_file_name = file.path(out_path,'/salmon_DESeq.feather')
write_feather(salmon_df, out_file_name)
message('Written: ', out_file_name)


# tximport kallisto abundance to gene count
abundance_table <- str_c(out_path,'/salmon_abundance.feather')
tximport(salmon_files_df$filename, 
        type = "salmon", 
        tx2gene = tx2gene, 
        reader = read_tsv,
        countsFromAbundance='lengthScaledTPM') %>%
    .$count %>%
    data.frame() %>%
    set_names(salmon_files_df$samplename) %>%
    rownames_to_column('id') %>%
    mutate(map_type = 'salmon') %>%
    write_feather(abundance_table)
message('Written: ', abundance_table)

