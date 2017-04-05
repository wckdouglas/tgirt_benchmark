#!/usr/bin/env Rscript

library(stringr)
library(readr)
library(dplyr)

bed_path <- '/stor/work/Lambowitz/ref/RNASeqConsortium'
gene_bed <- bed_path %>%
    str_c('/genes.bed') %>%
    read_tsv(col_names=c('chrom','start','end','gene_name','score','strand','bio_type','gene_id'),
             col_type = 'ciiccccc') %>%
    filter(!is.na(chrom)) %>%
    tbl_df

df <- read_tsv('/stor/work/Lambowitz/ref/human_transcriptome/transcripts.tsv') %>%
     right_join(gene_bed) %>%
     filter(!is.na(t_id) | bio_type=='tRNA')


filter_bed <- function(type_pattern, filename){
    bed_file_name <- str_c(bed_path, '/', filename,'.bed')
    if (!grepl('genes_no', filename)){
        gene_bed %>% 
            filter(grepl(type_pattern, bio_type)) %>%
            write_tsv(bed_file_name, col_names = F)
    }else{
        gene_bed %>%
            filter(!grepl(type_pattern, bio_type)) %>%
            write_tsv(bed_file_name, col_names=F)
    }
    message('Written ', bed_file_name)
}

patterns <- c('miRNA|misc_RNA|snoRNA|snRNA',
              'protein_coding',
             'miRNA|misc_RNA|snoRNA|snRNA|tRNA',
             '5S_rRNA|18S_rRNA|5.8S_rRNA|28S_rRNA',
             '18S_rRNA|28S_rRNA|5.8S_rRNA|5S_rRNA|miRNA|misc_RNA|rRNA|snoRNA|snRNA|tRNA',
             'rRNA',
             'miRNA|misc_RNA|rRNA|snoRNA|snRNA|tRNA')
filenames <- c('sncRNA_no_tRNA',
               'protein',
              'sncRNA_x_protein.bed',
              'rDNA', 
              'genes_no_sncRNA_rRNA_tRNA',
              'rRNA_for_bam_filter',
              'sncRNA_rRNA_for_bam_filter')
map2(patterns, filenames,filter_bed)
