#!/usr/bin/env Rscript

library(pROC)
library(readr)
library(dplyr)
library(feather)
library(cowplot)
library(stringr)
library(tidyr)
library(purrr)

project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'
# read gene table
figurepath <- str_c(project_path, '/figures')
ercc_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/ercc_annotation.tsv' %>%
    read_tsv()  %>%
    mutate(group = str_c(fold,1,sep=':')) %>%
    tbl_df

# read all tables ====================================================

df <- file.path(project_path, 'DEgenes') %>%
    list.files(path=., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|_[0-9]+|fc_|bias|align',.)] %>%
    map_df(read_feather) %>%
    mutate(id = str_replace(id, '_gene$','')) %>%
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
                grepl('Conventional',.$map_type) ~ "HISAT2+featureCounts",
                grepl('Customized', .$map_type) ~ "TGIRT-map",
                TRUE~ .$map_type)) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) %>%
    tbl_df

error_df <- df %>% 
    filter(!is.na(error)) %>%
    mutate(error = abs(error)) %>%
    group_by(map_type)  %>%
    summarize(me = mean(error),
              se = sd(error)) %>%
    ungroup() %>%
    mutate(y_coor = 1) %>%
    mutate(error_str = str_c(signif(me,2),' ± ',signif(se,2)))

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

all_comparison <- gtools::permutations(n=4,r=2,
                                       v=unique(as.character(rmse_df$map_type)),
                                       repeats.allowed=F)
get_all_p <-function(x1, x2, rmse_df){
    compare <- str_c(x1, x2, sep='|')
    p <- rmse_df %>% 
        filter(map_type == x1 | map_type == x2) %>%
        wilcox.test(rmse~map_type,paired=T, data=.) %>%
        .$p.value
    return(data.frame(p = p, comparison = compare))        
}
pval_df <- all_comparison %>%
    data.frame() %>%
    mutate(wilcox_p = map2(X1, X2, get_all_p, rmse_df)) %>%
    unnest(wilcox_p) %>%
    tbl_df

error_table <- str_c(figurepath, '/roc_rmse.csv')
rmse_df %>% 
    select(-y_coor) %>% 
    mutate(rmse = signif(rmse,3))%>% 
    spread(map_type, rmse) %>%
    rename(`Differentially-expressed group` = group) %>%
    write_csv(error_table) 

ercc_de_p<-ggplot()+
    geom_point(data=df %>% 
                inner_join(error_df) %>% 
                mutate(map_type = str_c(map_type,'\n|Error| = ',error_str)) %>%
                mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
                mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) , 
                        aes(x = log2(av_exp), 
                            y = log2FoldChange_AB, 
                            color = group)) +
    geom_hline(data=df %>% select(group, log2fold),
               aes(yintercept = log2fold, color = group)) +
#    geom_text(data=error_df, x= 4,
#              aes(label = str_c('Error: ',error), color = group, y = 5 - y_coor * 0.27)) +
    facet_grid(~map_type) +
    labs(x = 'Average expression (log2)', 
         y = 'Fold change between sample AB (log2)', 
         color = 'ERCC\nDE group') +
    scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"))

colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
error_plot <- ggplot(data=df, aes(fill = map_type, y=error, x = group)) +
    geom_violin()+
    labs(y =  expression(paste(Delta~'log2(fold change)')), 
         x = 'Expected ERCC differentially-expressed group (Sample B: Sample A)', 
         fill = ' ') +
    geom_hline(yintercept = 0, linetype=2, alpha=0.6) + 
    scale_fill_manual(values=colors)

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
    mutate(map_type = str_c(map_type, ' (AUC: ',auc,')')) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type))

roc_p <- ggplot(data=roc_df %>% arrange(tpr), aes(x = fpr, y = tpr, color = map_type)) +
    geom_line(size=2)+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.5,0.25)) +
    geom_abline(slope = 1, intercept = 0, linetype=2)+
    scale_color_manual(values=colors)


p <- plot_grid(error_plot, roc_p, 
               ncol=1, labels = letters[1:2], 
               label_size=20,
               align='v', axis='l')
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/ercc_roc.pdf')
save_plot(p, file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)

ercc_roc_table <- str_c(figurepath, '/ercc_roc.csv')
roc_df %>% 
    mutate(map_type = str_replace(map_type,' \\(AUC: 0.[0-9]+\\)','')) %>%
    select(map_type, auc) %>%
    dplyr::rename(auc_ercc = auc) %>%
    distinct() %>%
    write_csv(ercc_roc_table)
