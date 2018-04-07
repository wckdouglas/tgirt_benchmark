#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(purrr)
library(feather)


get_mix <- function(samplename){
    ifelse(grepl('^ref',samplename), str_sub(samplename,4,4), str_sub(samplename,8,8))
}

get_sample_number <- function(samplename){
    id <- ifelse(grepl('^ref',samplename), str_sub(samplename,5,5), str_sub(samplename,10,10))
    return(as.numeric(id))
}

gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv'  %>%
    read_tsv() %>%
    dplyr::select(gene_id, name,type) %>% 
    unique %>%
    dplyr::rename(id = gene_id)

ercc_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/ercc_annotation.tsv' %>%
    read_tsv() %>%
    inner_join(gene_file)



#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
    .[!grepl('_[0-9]+',.)] %>%
    map_df(read_feather) %>%
    gather(samplename, abundance, -id, -map_type) %>%
    inner_join(gene_file) %>%
    filter(type=='ERCC') %>%
    mutate(sample_mix = get_mix(samplename)) %>%
    mutate(sample_id = get_sample_number(samplename)) %>%
    inner_join(ercc_file) %>%
    mutate(map_type = case_when(
                                grepl('conventional',.$map_type) ~ "HISAT2+featureCounts",
                                grepl('customized', .$map_type) ~ "TGIRT-map",
                                TRUE~ .$map_type)) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type))

df %>% 
    mutate(samplename = str_replace(samplename,'_[123]','')) %>%
    group_by(map_type, id, samplename) %>% 
    summarize(abundance=mean(abundance)) %>% 
    ungroup %>% 
    filter(abundance>0)%>% 
    group_by(samplename, map_type) %>% 
    summarize(number_og_genes=n())

lm_df <- df %>%
    filter(group %in% c('A','B')) %>%
    mutate(conc = ifelse(sample_mix=='A',mix1, mix2)) %>%
    tbl_df

formula = y~x
ercc_lm <- ggplot(data = lm_df, aes(x = log2(conc), y = log2(abundance))) +
    geom_point(alpha=0.6) +
    geom_smooth(method='lm', formula = formula) +
    ggpmisc::stat_poly_eq(formula = formula, parse = TRUE) +
    facet_grid(group~map_type)+
    panel_border() +
    scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")) +
    labs(x = 'Concentration [log2(amol/ul)]', y = 'TPM (log2)') +
    theme(legend.position = 'none') 

r2_df <- lm_df %>%
    mutate(log2_abundance = log2(abundance)) %>%
    filter(!is.infinite(log2_abundance)) %>%
    mutate(log2_conc = log2(conc)) %>%
    group_by(map_type, group, sample_id) %>%
    nest() %>%
    mutate(model = map(data, ~lm(log2_abundance~log2_conc, data=.))) %>%
    mutate(r2 = map(model, broom::glance)) %>%
    unnest(r2) %>% 
    tbl_df
pval <- r2_df %>% 
    mutate(map_type = factor(map_type)) %>% 
    kruskal.test(r.squared~map_type, data=.) %>%
    .$p.value

ercc_r2 <- ggplot(data=r2_df, aes(x = map_type, color=map_type, y = r.squared, shape=group)) +
    geom_jitter(size=3, width = 0.1) +
    labs(x = ' ', y = expression(R^2), parse=T, shape = 'Sample mix', color = 'Pipeline') +
    theme(axis.title.y = element_text(angle=0)) +
    ylim(0.9,1) +
    scale_colour_manual(values=RColorBrewer::brewer.pal(8, "Dark2"))

p <- plot_grid(ercc_lm, ercc_r2, ncol=1, labels = letters[1:2])
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/ercc_figure.png')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)


