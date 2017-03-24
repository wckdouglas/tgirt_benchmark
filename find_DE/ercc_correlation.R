#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)

filter_func <- function(d, label){
    d %>%
    filter(grepl('ERCC',id))  %>%
    gather(samplename, abundance, -id) %>%
    mutate(prep = label) %>%
    tbl_df 
}

get_mix <- function(samplename){
    ifelse(grepl('^ref',samplename), str_sub(samplename,4,4), str_sub(samplename,8,8))
}

get_sample_number <- function(samplename){
    id <- ifelse(grepl('^ref',samplename), str_sub(samplename,5,5), str_sub(samplename,10,10))
    return(as.numeric(id))
}


ercc_file <- '/stor/work/Lambowitz/ref/RNASeqConsortium/ercc/ercc_table.tsv' %>%
    read_tsv()  %>%
    tbl_df

project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking'
kallisto_count <- project_path %>%
    str_c('/alignment_free/countFiles/alignment_free.tsv')  %>%
    read_tsv() %>%
    dplyr::rename(id = target_id)  %>%
    filter_func('Alignment Free') %>%
    tbl_df
    
multimap_count <- project_path %>%
    str_c('/genome_mapping/pipeline7_counts/RAW/combined_gene_count.tsv')  %>%
    read_tsv() %>%
    filter_func('W/ mutlimap') %>%
    tbl_df

    
genome_map_count <- str_c(project_path, '/genome_mapping/old_pipeline/countsData.tsv')   %>%
    read_tsv() %>%
    select(grep('id|ref',names(.))) %>%
    filter_func('W/o mutlimap') %>%
    tbl_df



df <- purrr::reduce(list(kallisto_count, multimap_count, genome_map_count), rbind) %>%
    mutate(sample_mix = get_mix(samplename)) %>%
    mutate(sample_id = get_sample_number(samplename)) %>%
    inner_join(ercc_file)

lm_df <- df %>%
    filter(group %in% c('A','B')) %>%
    mutate(conc = ifelse(sample_mix=='A',mix1, mix2)) %>%
    tbl_df

formula = y~x
ercc_lm <- ggplot(data = lm_df, aes(x = log2(conc), y = log2(abundance))) +
    geom_point(alpha=0.6) +
    geom_smooth(method='lm', formula = formula) +
    ggpmisc::stat_poly_eq(formula = formula, parse = TRUE) +
    facet_grid(group~prep)+
    panel_border()

r2_df <- lm_df %>%
    mutate(log2_abundance = log2(abundance)) %>%
    filter(!is.infinite(log2_abundance)) %>%
    mutate(log2_conc = log2(conc)) %>%
    group_by(prep, group, sample_id) %>%
    nest() %>%
    mutate(model = map(data, ~lm(log2_abundance~log2_conc, data=.))) %>%
    mutate(r2 = map(model, broom::glance)) %>%
    unnest(r2) %>% 
    tbl_df

ercc_r2 <- ggplot(data=r2_df, aes(x = prep, color=prep, y = r.squared, shape=group)) +
    geom_jitter(size=3, width = 0.1) +
    labs(x = ' ', y = expression(R^2), parse=T, shape = 'Sample mix', color = 'Pipeline') +
    theme(axis.title.y = element_text(angle=0)) +
    ylim(0.9,1)

p <- plot_grid(ercc_lm, ercc_r2, ncol=1)
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/ercc_figure.png')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)


