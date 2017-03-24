#!/usr/bin/evn Rscript

library(feather)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(purrr)

taqman <- '/stor/scratch/Lambowitz/cdw2854/bench_marking/maqc/taqman_fc_table.feather' %>%
    read_feather() %>%
    rename(id = ensembl_gene_id) %>%
    mutate(real_FC = ifelse(abs(logFC_AB)>0.5, 'DE','notDE')) %>%
    rename(taqman_fc_AB = logFC_AB)

project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking'
alignment_free_df <- project_path %>%
    str_c('/alignment_free/countFiles/sleuth_results.feather') %>%
    read_feather() %>%
    mutate(logFC = -b) %>%
    dplyr::rename(id = gene_id) %>%
    select(id,sample_base,sample_test, name,type,qval, mean_obs, logFC, pval) %>%
    gather(variable, value, -id,-name,-type,-sample_base,-sample_test) %>%
    mutate(variable = case_when(
        .$variable == 'qval' ~ 'adj.P.Val',
        .$variable == 'mean_obs' ~ 'AveExpr',
        .$variable == 'logFC' ~ 'logFC',
        .$variable == 'pval' ~ 'P.Val'
    ))  %>%
    mutate(variable = str_c(variable,'_', sample_base, sample_test)) %>%
    select(-sample_base,-sample_test) %>%
    spread(variable, value)  %>%
    select(-name,-type) %>%
    mutate(map_type = 'Kallisto') %>%
    tbl_df

genome_df <- project_path %>%
    str_c('/genome_mapping/pipeline7_counts/deseq_genome.feather') %>%
    read_feather() %>%
    tbl_df
df <- rbind(genome_df, alignment_free_df) %>%
    inner_join(taqman)

## plotting roc curve
roc_data <- function(p_cut, de_df){
    de_df %>%
        mutate(DE = ifelse(adj.P.Val_AB<p_cut, 'DE','notDE')) %>%
        group_by(map_type,DE, real_FC) %>%
        summarize(count = n())  %>%
        ungroup() %>%
        mutate(p = p_cut) 
}
# 23 non DE in ERCC, 69 DE
roc_df <- df %>% 
    mutate(adj.P.Val_AB = ifelse(is.na(adj.P.Val_AB),1,adj.P.Val_AB)) %>%
    map_df(seq(0,1,0.001), roc_data,.)  %>%
    mutate(label = str_c(DE,'-',real_FC)) %>%
    mutate(label = case_when(
        .$label == 'notDE-DE' ~ 'FN',
        .$label == 'notDE-notDE' ~ 'TN',
        .$label == 'DE-DE' ~ 'TP',
        .$label == 'DE-notDE' ~ 'FP'
    )) %>%
    select(-DE, -real_FC) %>%
    spread(label,count, fill=0)%>%
    mutate(TPR = TP/(TP+FN)) %>%
    mutate(FPR = FP/(TN+FP)) %>%
    mutate(slope = TPR/FPR) %>%
    tbl_df

label_df <-  roc_df %>% 
    group_by(map_type) %>%
    top_n(1,slope) %>%
    do(head(.,1)) %>%
    ungroup() %>%
    mutate(label = str_c('p = ',signif(p,3))) %>%
    tbl_df

gene_roc_p <- ggplot(data=roc_df, aes(x = FPR, y = TPR, color = map_type)) +
    geom_point() +
    geom_line()+
    ggrepel::geom_label_repel(data=label_df, 
                              aes(x = FPR, y = TPR, label = label, color = map_type))+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.75,0.25)) +
    geom_abline(slope = 1, intercept = 0)


rmsd_df <- df %>%
    select(map_type, logFC_AB, taqman_fc_AB) %>%
    mutate(sq_error = (logFC_AB - taqman_fc_AB)^2) %>%
    group_by(map_type) %>%
    summarize(rmse = sqrt(mean(sq_error, na.rm=T))) %>%
    ungroup() 

rmse_bar <- ggplot(data=rmsd_df, aes(x = map_type, y = rmse, fill=map_type)) +
    geom_bar(stat='identity') +
    labs(x = ' ', y = 'RMSE (RNA-seq vs TaqMan)') +
    theme(legend.position = 'none')

p <- plot_grid(rmse_bar, gene_roc_p)
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/taqman_roc.png')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
    
