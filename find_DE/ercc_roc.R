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
ercc_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/ercc_annotation.tsv' %>%
    read_tsv()  %>%
    mutate(group = str_c(fold,1,sep=':')) %>%
    tbl_df

# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'

df <- file.path(project_path, 'DEgenes') %>%
    list.files(path=., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|_[0-9]+',.)] %>%
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
    mutate(map_type = case_when(
                grepl('Conventional',.$map_type) ~ "HISAT2+FeatureCounts",
                grepl('Customized', .$map_type) ~ "TGIRT-map",
                TRUE~ .$map_type)) %>%
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
pval <- rmse_df %>% 
    mutate(map_type = factor(map_type))%>% 
    friedman.test(rmse~map_type|group,data=.) %>%
    .$p.value
sal_vs_kall <- rmse_df %>% filter(grepl('Ka|Sa',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)
sal_vs_conv <- rmse_df %>% filter(grepl('HI|Sa',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)
ka_vs_conv <- rmse_df %>% filter(grepl('HI|Ka',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)
ka_vs_cust <-rmse_df %>% filter(grepl('TG|Ka',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)
sal_vs_cust <-rmse_df %>% filter(grepl('TG|Sa',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)
con_vs_cust <-rmse_df %>% filter(grepl('TG|HI',map_type))%>% wilcox.test(rmse~map_type,paired=T, data=.)

ercc_de_p<-ggplot()+
    geom_point(data=df, aes(x = log2(av_exp), 
                            y = log2FoldChange_AB, 
                            color = group)) +
    geom_hline(data=df %>% select(group, log2fold),
               aes(yintercept = log2fold, color = group)) +
    geom_text(data=rmse_df, x= 10,
              aes(label = str_c('RMSE: ',signif(rmse,3)), color = group, y = 3.5 - y_coor * 0.25)) +
    facet_grid(~map_type) +
    labs(x = 'Average expression (log2)', 
         y = 'Fold change between sample AB (log2)', 
         color = 'ERCC\nDE group') +
    scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"))

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


p <- plot_grid(ercc_de_p, roc_p, ncol=1, labels = letters[1:2], label_size=20)
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/ercc_roc.pdf')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)

ercc_roc_table <- str_c(figurepath, '/ercc_roc.csv')
roc_df %>% 
    mutate(map_type = str_replace(map_type,' \\(AUC: 0.[0-9]+\\)','')) %>%
    select(map_type, auc) %>%
    rename(auc_ercc = auc) %>%
    distinct() %>%
    write_csv(ercc_roc_table)