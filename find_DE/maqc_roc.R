#!/usr/bin/evn Rscript

library(feather)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(AUC)
library(purrr)

taqman <- '/stor/scratch/Lambowitz/cdw2854/bench_marking/maqc/taqman_fc_table.feather' %>%
    read_feather() %>%
    dplyr::rename(id = ensembl_gene_id) %>%
    mutate(real_FC = ifelse(abs(logFC_AB)>0.5, 'DE','notDE')) %>%
    rename(taqman_fc_AB = logFC_AB)

# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance',.)] %>%
    map_df(read_feather) %>%
    gather(variable, value, -id, -map_type, - comparison) %>%
    filter(grepl('pvalue|padj|baseMean|log2FoldChange', variable)) %>%
    mutate(sample1 = str_sub(comparison, 1, 1)) %>%
    mutate(sample2 = str_sub(comparison, 6, 6)) %>%
    mutate(variable = str_c(variable,'_', sample1, sample2)) %>%
    select(-sample1, -sample2,- comparison) %>%
    spread(variable, value) %>%
    inner_join(taqman)

## plotting roc curve
roc_data <- function(p_cut, de_df){
    de_df %>%
        mutate(DE = ifelse(pvalue_AB<p_cut, 'DE','notDE')) %>%
        group_by(map_type,DE, real_FC) %>%
        summarize(count = n())  %>%
        ungroup() %>%
        mutate(p = p_cut) 
}
# 23 non DE in ERCC, 69 DE
roc_df <- df %>% 
    mutate(pvalue_AB = ifelse(is.na(pvalue_AB),1,pvalue_AB)) %>%
    mutate(label = factor(ifelse(real_FC=='DE',0,1))) %>%
    group_by(map_type) %>%
    nest() %>%
    mutate(roc_model = map(data, ~AUC::roc(.$padj_AB, .$label))) %>%
    mutate(tpr = map(roc_model, ~.$tpr)) %>%
    mutate(fpr = map(roc_model, ~.$fpr))
    
AUC_df <- roc_df %>%
    mutate(auc = map_dbl(roc_model, AUC::auc))%>%
    select(auc, map_type) %>%
    mutate(auc = signif(auc,3)) 
    

roc_plot_df <- roc_df %>%
    unnest(tpr, fpr) %>%
    inner_join(AUC_df) %>%
    mutate(map_type = str_c(map_type, ' (AUC: ',auc,')'))

    

gene_roc_p <- ggplot(data=roc_plot_df, aes(x = fpr, y = tpr, color = map_type)) +
    geom_line()+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.75,0.25)) +
    geom_abline(slope = 1, intercept = 0)


rmsd_df <- df %>%
    select(map_type, log2FoldChange_AB, taqman_fc_AB) %>%
    mutate(sq_error = (log2FoldChange_AB - taqman_fc_AB)^2) %>%
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
message('Written: ', figurename)
    
