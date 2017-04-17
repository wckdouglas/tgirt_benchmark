#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(tidyr)
library(feather)



#======set up change type funciton
ncRNA=c("sense_intronic","3prime_overlapping_ncRNA",'processed_transcript',
        'sense_overlapping','Other_lncRNA', 'macro_lncRNA','non_coding',
        'lincRNA','bidirectional_promoter_lncRNA', 'ribozyme')
smncRNA=c('misc_RNA','snRNA','piRNA','scaRNA','sRNA','scRNA')
large_rRNA=c('28S_rRNA','18S_rRNA')
small_rRNA=c('rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA')
protein_coding = c('protein_coding','TR','IG')

changeType <- function(x){
    if (x %in% ncRNA){
        'Other ncRNA'
    }else if (grepl('TR|IG|protein|pseudogene',x)){
        'Protein coding'
    }else if (grepl('Mt_',x)){
        'Mt'
    }else if (grepl('tRNA',x)){
        'tRNA'
    }else if (x %in% small_rRNA){
        '5/5.8S rRNA'
    }else if (x %in% large_rRNA){
        '18/28S rRNA'
    }else if (x %in% smncRNA){
        'Other sncRNA'
    }else if (x =='antisense'){
        'Antisense'
    }else {
        x
    }
}


tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}


#read files
get_mix <- function(samplename){
    ifelse(grepl('^ref',samplename), str_sub(samplename,4,4), str_sub(samplename,8,8))
}

get_sample_number <- function(samplename){
    id <- ifelse(grepl('^ref',samplename), str_sub(samplename,5,5), str_sub(samplename,10,10))
    return(as.numeric(id))
}

gene_file <- '/stor/work/Lambowitz/ref/human_transcriptome/transcripts.tsv'  %>%
    read_tsv() %>%
    select(gene_id, name,type) %>% 
    unique %>%
    dplyr::rename(id = gene_id) %>%
    mutate(type = sapply(type,changeType)) %>%
    tbl_df

gene_length_df <- '/stor/work/Lambowitz/ref/human_transcriptome/genes.length' %>%
    read_tsv
 

#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
alignment_free <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
    map_df(read_feather) 

# read genome mapping count files
files <- c(file.path(project_path, 'genome_mapping/conventional/feature_counts.tsv'),
           file.path(project_path,'genome_mapping/pipeline7/Counts/RAW/combined_gene_count.tsv'))
labels <- c('conventional','customized')
genome_df <- map2(files, labels, function(x,y) read_tsv(x) %>% 
                      mutate(map_type=y) %>%
                      set_names(str_replace_all(names(.),'-','_'))) %>%
    purrr::reduce(rbind) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]+$','') ,id)) %>%
    tbl_df

df <- rbind(genome_df, alignment_free) %>%
    gather(samplename, abundance, -id, -map_type) %>%
    inner_join(gene_file) %>%
    group_by(type, map_type, samplename) %>%
    summarize(read_count = sum(abundance)) %>%
    ungroup() %>%
    group_by(map_type, samplename )%>%
    do(data_frame(
        percentage_count = .$read_count/sum(.$read_count) * 100,
        type=.$type
    )) %>%
    ungroup %>%
    mutate(mix = get_mix(samplename)) %>%
    mutate(replicate = get_sample_number(samplename)) %>%
    mutate(x_name = str_c(mix, replicate)) %>%
    mutate(type = forcats::fct_reorder(type, percentage_count, sum)) %>%
    tbl_df


#make color
library(RColorBrewer)
n <- length(unique(df$type))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(33)

dist_p <-ggplot(data = df, aes(x=x_name, y = percentage_count, fill=type)) + 
    geom_bar(stat='identity') + 
    facet_grid(map_type~mix, scale='free_x') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12,'Set3'))
    