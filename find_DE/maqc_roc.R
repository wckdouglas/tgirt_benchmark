#!/usr/bin/evn Rscript

library(pROC)
library(feather)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(purrr)
library(viridis)

taqman <- '/stor/work/Lambowitz/cdw2854/bench_marking/maqc/taqman_fc_table.feather' %>%
    read_feather() %>%
    dplyr::rename(id = ensembl_gene_id) %>%
    mutate(real_FC_AB = ifelse(abs(logFC_AB)>0.5, 'DE','notDE')) %>%
    mutate(real_FC_CD = ifelse(abs(logFC_CD)>0.5, 'DE','notDE')) %>%
    rename(taqman_fc_AB = logFC_AB) %>%
    rename(taqman_fc_CD = logFC_CD) %>%
    dplyr::select(id, real_FC_AB, real_FC_CD, taqman_fc_AB, taqman_fc_CD) %>%
    gather(variable, value, -id) %>% 
    mutate(comparison = str_sub(variable,-2,-1)) %>% 
    mutate(variable = str_sub(variable,1,-4))    %>%
    spread(variable, value) %>%
    mutate(taqman_fc = as.numeric(taqman_fc))

# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm',.)] %>%
    map_df(read_feather) %>%
    gather(variable, value, -id, -map_type, - comparison) %>%
    filter(grepl('pvalue|log2FoldChange', variable)) %>%
    mutate(comparison = str_replace(comparison,' vs ','')) %>%
    spread(variable, value) %>%
    inner_join(taqman) %>%
    tbl_df
    
df <- df %>% 
    filter(!is.na(log2FoldChange)) %>%
    group_by(id) %>% 
    summarize(occur = n()) %>% 
    ungroup %>% 
    filter(occur==8) %>% 
    inner_join(df) %>%
    select(-occur)

# 23 non DE in ERCC, 69 DE
roc_df <- df %>% 
    mutate(pvalue = ifelse(is.na(pvalue),1,pvalue)) %>%
    group_by(map_type) %>%
    nest() %>%
    mutate(roc_model = map(data, ~pROC::roc(.$real_FC, .$pvalue))) %>%
    mutate(tpr = map(roc_model, ~.$sensitivities)) %>%
    mutate(fpr = map(roc_model, ~1 - .$specificities))
    
AUC_df <- roc_df %>%
    mutate(auc = map_dbl(roc_model, pROC::auc))%>%
    dplyr::select(auc, map_type) %>%
    mutate(auc = signif(auc,3)) 
    

roc_plot_df <- roc_df %>%
    unnest(tpr, fpr) %>%
    inner_join(AUC_df) %>%
    mutate(map_type = str_c(map_type, ' (AUC: ',auc,')')) %>%
    arrange(tpr)

gene_roc_p <- ggplot(data=roc_plot_df, aes(x = fpr, y = tpr, color = map_type)) +
    geom_line()+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.7,0.25)) +
    geom_abline(slope = 1, intercept = 0)


rmsd_df <- df %>%
    dplyr::select(map_type, log2FoldChange, taqman_fc) %>%
#    mutate(sq_error = sqrt((log2FoldChange - taqman_fc)^2)) %>%
    mutate(sq_error = log2FoldChange - taqman_fc) %>%
    ungroup() 
pval <- kruskal.test(sq_error~factor(map_type), data=rmsd_df)$p.value

rmse_bar <- ggplot(data=rmsd_df, aes(x = map_type, y = abs(sq_error), fill = map_type)) +
    geom_violin() +
    labs(x = ' ', y = expression(paste(Delta~'log2(fold change)'))) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_blank())  

p <- plot_grid(rmse_bar, gene_roc_p, labels = letters[1:2], label_size=20)
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/taqman_roc.pdf')
save_plot(p, file=figurename,  base_width=12, base_height=7) 
message('Written: ', figurename)
    
missing_gene <- df %>% 
#    filter(comparison == 'AB') %>% 
    filter(comparison == 'CD') %>% 
    select(-pvalue, - real_FC, - log2FoldChange) %>% 
    spread(map_type, taqman_fc) %>% 
    filter(is.na(Conventional_pipeline))

df %>% group_by(comparison, map_type) %>% summarize(n())
