#!/usr/bin/env Rscript

# download file:   wget ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4097/soft/GPL4097_family.soft.gz
# zcat GPL4097_family.soft.gz | grep -v '^!\|^#\|^\^' | awk 'NF==5' > GPL4097_annot.tsv


library(stringr)
library(biomaRt)
library(purrr)
library(readr)
library(feather)
library(dplyr)
library(tidyr)
library(GEOquery)

work_path <- "/stor/work/Lambowitz/cdw2854/bench_marking/maqc"
setwd(work_path)

# download id conversion table
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- c('refseq_mrna','ensembl_gene_id')
id_conversions <- getBM(attributes=attributes, 
      mart = ensembl) 
id_file <- str_c(work_path,'/refseq2ensembl.feather')
write_feather(id_conversions, id_file)
message('Written ', id_file)

# download taqman results
gset <- getGEO("GSE5350", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL4097", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
probe_table <- read_tsv('GPL4097_annot.tsv') %>% 
    mutate(GB_LIST=sapply(GB_LIST, function(x) str_split(x,',')[[1]][1]))


column_label <- gset@phenoData@data$title
maqc_expr_data <- gset@assayData$exprs %>%
    tbl_df %>%
    set_names(as.character(column_label)) %>%
    dplyr::select(grep('[ABCD][0-9]$',names(.))) %>%
    mutate(refseq_mrna = probe_table$GB_LIST) %>%
    gather(exp_sample, expr_value,-refseq_mrna) %>%
    mutate(sample_mix = str_sub(exp_sample,12,12)) %>%
    dplyr::select(-exp_sample) %>%
    group_by(sample_mix, refseq_mrna) %>%
    summarize(expr_value  = mean(expr_value)) %>%
    ungroup() %>%
    inner_join(id_conversions) %>%
    spread(sample_mix, expr_value) %>%
    mutate(logFC_AB=log2(B/A)) %>%
    mutate(logFC_CD=log2(D/C))
expr_data <- str_c(work_path, '/taqman_fc_table.feather')
write_feather(maqc_expr_data,expr_data)
message('Written: ', expr_data)


