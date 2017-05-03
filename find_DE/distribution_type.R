#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(tidyr)
library(feather)
library(purrr)


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


count_to_tpm <- function(counts, lengths) {
    rate <- counts / lengths
    out_tpm <- rate / sum(rate,na.rm=T) * 1e6
    return(out_tpm)
}


#read files
get_mix <- function(samplename){
    ifelse(grepl('^ref',samplename), str_sub(samplename,4,4), str_sub(samplename,8,8))
}

get_sample_number <- function(samplename){
    id <- ifelse(grepl('^ref',samplename), str_sub(samplename,5,5), str_sub(samplename,10,10))
    return(as.numeric(id))
}

gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv'  %>%
    read_tsv() %>%
    select(gene_id, name,type) %>% 
    unique %>%
    dplyr::rename(id = gene_id) %>%
    mutate(type = sapply(type,changeType)) %>%
    mutate(type= ifelse(grepl('MT-T',name),'tRNA',type)) %>% 
    mutate(id = ifelse(type=='tRNA', name, id)) %>%
    tbl_df

#gene_length_df <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/genes.length' %>%
#    read_tsv  %>%
#    unique()
 

#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
alignment_free <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
    map_df(read_feather)  %>%
    gather(samplename, abundance, -id, -map_type) %>%
    tbl_df

# read genome mapping count files
files <- c(file.path(project_path, 'genome_mapping/Trim/conventional/counts/feature_counts.tsv'),
           file.path(project_path,'genome_mapping/Counts/RAW/combined_gene_count.tsv'))
labels <- c('conventional','customized')
genome_df <- map2(files, labels, function(x,y) read_tsv(x) %>% 
                      mutate(map_type=y) %>%
                      set_names(str_replace_all(names(.),'-','_'))) %>%
    purrr::reduce(rbind) %>%
#    inner_join(gene_length_df) %>%
    mutate(id = str_replace(id,'\\-$',''))%>%
    gather(samplename, abundance, -id,-map_type)%>%#, -gene_length) %>%
    mutate(id = ifelse(!grepl('ERCC',id),str_replace(id, '\\-[0-9]+$',''),id)) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]+$','') ,id)) %>%
    group_by(id, samplename, map_type) %>%
    summarize(
        abundance = sum(abundance)#, 
#       gene_length = mean(gene_length)
    ) %>%
    ungroup() %>%
#    mutate(tpm = count_to_tpm(abundance, gene_length)) %>%
#    select(-abundance, - gene_length) %>%
#    dplyr::rename(abundance=tpm) %>%
    tbl_df

merge_df <- rbind(genome_df, alignment_free) %>%
    inner_join(gene_file) %>%
    tbl_df

merge_df <-  merge_df %>% 
    filter(samplename == 'Sample_A_1') %>%
    group_by(id,samplename) %>% 
    summarize(count = n(), 
              types = str_c(map_type,collapse=',')) %>% 
    filter(count == 4) %>%
    select(-samplename, -types, count) %>%
    inner_join(merge_df) %>%
    ungroup() 

#make color
gene_count_p <- ggplot(data = merge_df %>% 
                filter(abundance > 0) %>%
                group_by(type,samplename, map_type) %>%
                summarize(gene_count = n()) %>%
                ungroup() %>%
                mutate(samplename = str_replace(samplename, 'Sample_','')) %>%
                mutate(samplename = str_replace(samplename, '_','')) %>%
                mutate(type = forcats::fct_reorder(type, gene_count , sum)),
        aes(x = samplename, y = gene_count, fill = type)) +
    geom_bar(stat='identity') +
    facet_grid(~map_type) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12,'Set3')) +
    labs(x = ' ', y = 'Number of detected genes', fill= ' ')


    
df <- merge_df %>%
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
    mutate(type = forcats::fct_reorder(type, percentage_count , sum))%>%
    tbl_df


dist_p <-ggplot(data = df, aes(x=x_name, y = percentage_count, fill=type)) + 
    geom_bar(stat='identity') + 
    facet_grid(~map_type) +
    labs(x = ' ', y = '% count', fill= ' ') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12,'Set3'))

p<-plot_grid(gene_count_p, dist_p, 
             labels=letters[1:2], ncol=1)
    

plot_pairs <- function(rna_type){
    p <- merge_df %>% 
        filter(grepl('Sample_A',samplename)) %>%
        filter(grepl(rna_type,type)) %>% 
        group_by(map_type, id) %>%
        summarize(abundance = log2(mean(abundance)+1)) %>%
        ungroup()  %>%
        spread(map_type, abundance) %>%
        select(-id) %>% 
        GGally::ggpairs()
}

ps <- lapply(c('tRNA','snoRNA','Protein'),plot_pairs)
