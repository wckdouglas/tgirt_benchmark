#!/usr/bin/env Rscript

library(biomaRt)
library(dplyr)
library(readr)
library(purrr)

args<-commandArgs(TRUE)
out_bed_name <- args[1]
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bed_fields <- c('chromosome_name',
                'start_position',
                'end_position',
                'external_gene_name',
                'strand',
                'gene_biotype',
                'ensembl_gene_id')
getBM(attributes = bed_fields,
      mart = ensembl) %>%
    set_names(c('chrom','start','end','name','strand','biotype','id')) %>%
    mutate(score=0) %>%
    dplyr::select(chrom,start,end,name,score,strand, biotype, id) %>%
    arrange(chrom,start) %>%
    mutate(strand = ifelse(strand < 1, '-','+')) %>%
    filter(biotype!='TEC') %>%
    write_tsv(out_bed_name, col_names=F)
message('Downloaded ', out_bed_name)
    