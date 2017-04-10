#!/usr/bin/env Rscript

library(pROC)
library(readr)
library(dplyr)
library(feather)
library(cowplot)
library(stringr)
library(tidyr)
library(purrr)

# read gene table
ercc_file <- '/stor/work/Lambowitz/ref/RNASeqConsortium/ercc/ercc_table.tsv' %>%
    read_tsv()  %>%
    tbl_df

# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'


df <- file.path(project_path, 'DEgenes') %>%
    list.files(path=., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance',.)] %>%
    map_df(read_feather) %>%
    gather(variable, value, -id, -map_type, - comparison) %>%
    filter(grepl('ERCC',id)) %>%
    filter(grepl('pvalue|padj|baseMean|log2FoldChange', variable)) %>%
    mutate(sample1 = str_sub(comparison, 1, 1)) %>%
    mutate(sample2 = str_sub(comparison, 6, 6)) %>%
    mutate(variable = str_c(variable,'_', sample1, sample2)) %>%
    select(-sample1, -sample2,- comparison) %>%
    spread(variable, value) %>%
    inner_join(ercc_file,by='id') %>%
    select(map_type, id, log2FoldChange_AB, pvalue_AB, group,fold, mix1,mix2, padj_AB,label) %>%
    mutate(log2fold = log2(fold)) %>%
    mutate(av_exp = (mix1+mix2)/2) %>%
    mutate(padj_AB = ifelse(is.na(padj_AB),1,padj_AB)) %>%
    mutate(pvalue_AB = ifelse(is.na(pvalue_AB),1,pvalue_AB)) %>%
    mutate(error = log2FoldChange_AB-log2fold)  %>% 
    tbl_df

rmse_df <- df %>% 
    filter(!is.na(error)) %>%
    group_by(map_type, group)  %>%
    summarize(rmse = sqrt(mean(error^2))) %>%
    ungroup() %>%
    mutate(y_coor = 1) %>%
    group_by(map_type) %>%
    do(data_frame(
        rmse =.$rmse,
        group = .$group,
        y_coor = cumsum(.$y_coor)
    )) %>%
    ungroup()

ercc_de_p<-ggplot()+
    geom_point(data=df, aes(x = log2(av_exp), 
                            y = log2FoldChange_AB, 
                            color = group)) +
    geom_hline(data=df %>% select(group, log2fold),
               aes(yintercept = log2fold, color = group)) +
    geom_text(data=rmse_df, x= 10,
              aes(label = str_c('RMSE: ',signif(rmse,3)), color = group, y = 3.5 - y_coor * 0.25)) +
    facet_grid(~map_type) +
    labs(x = 'log2(average expression)', 
         y = 'log2(Fold change between AB)', 
         color = ' ')

# 23 non DE in ERCC, 69 DE
ercc_roc_df <- df %>% 
    mutate(padj_AB = ifelse(is.na(padj_AB),1,padj_AB)) %>%
    mutate(pvalue_AB = ifelse(is.na(pvalue_AB),1,pvalue_AB)) %>%
    group_by(map_type) %>%
    nest() %>%
    mutate(roc_model = map(data, ~pROC::roc(.$label,.$padj_AB))) %>%
    mutate(tpr = map(roc_model, ~.$sensitivities)) %>%
    mutate(fpr = map(roc_model, ~1 - .$specificities))

auc_df <- ercc_roc_df %>%
    mutate(auc = map_dbl(roc_model, pROC::auc)) %>%
    select(map_type, auc) %>%
    mutate(auc = signif(auc,3))

roc_df <- ercc_roc_df %>%
    unnest(tpr, fpr) %>%
    inner_join(auc_df) %>%
    mutate(map_type = str_c(map_type, ' (AUC: ',auc,')'))

roc_p <- ggplot(data=roc_df %>% arrange(tpr), aes(x = fpr, y = tpr, color = map_type)) +
    geom_line()+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.75,0.25)) +
    geom_abline(slope = 1, intercept = 0)


p <- plot_grid(ercc_de_p, roc_p, ncol=1)
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/ercc_roc.png')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)
