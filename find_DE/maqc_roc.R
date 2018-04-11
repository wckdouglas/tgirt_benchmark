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


project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'

taqman <- file.path(project_path,'maqc/taqman_fc_table.feather') %>%
    read_feather() %>%
    dplyr::rename(id = ensembl_gene_id) %>%
    mutate(real_FC_AB = ifelse(abs(logFC_AB)>0.5, 'DE','notDE')) %>%
    mutate(real_FC_CD = ifelse(abs(logFC_CD)>0.5, 'DE','notDE')) %>%
    dplyr::rename(taqman_fc_AB = logFC_AB) %>%
    dplyr::rename(taqman_fc_CD = logFC_CD) %>%
    dplyr::select(id, real_FC_AB, real_FC_CD, taqman_fc_AB, taqman_fc_CD) %>%
    gather(variable, value, -id) %>% 
    mutate(comparison = str_sub(variable,-2,-1)) %>% 
    mutate(variable = str_sub(variable,1,-4))    %>%
    spread(variable, value) %>%
    mutate(taqman_fc = as.numeric(taqman_fc))

# read all tables ====================================================
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|_[0-9]+|bias|aligned',.)] %>%
    map_df(read_feather) %>%
    gather(variable, value, -id, -map_type, - comparison) %>%
    filter(grepl('pvalue|log2FoldChange', variable)) %>%
    mutate(comparison = str_replace(comparison,' vs ','')) %>%
    spread(variable, value) %>%
    inner_join(taqman) %>%
    mutate(map_type = case_when(
                grepl('Conventional',.$map_type) ~ "HISAT2+featureCounts",
                grepl('Customized', .$map_type) ~ "TGIRT-map",
                TRUE~ .$map_type)) %>%
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
    arrange(tpr) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) 

colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
gene_roc_p <- ggplot(data=roc_plot_df, aes(x = fpr, y = tpr, color = map_type)) +
    geom_line(size=2, alpha=0.8)+
    labs(x = 'False positive rate', y = 'True positive rate', color = ' ') +
    theme(legend.position = c(0.6,0.25)) +
    geom_abline(slope = 1, intercept = 0, linetype=2) +
    scale_color_manual(values=colors)


rmsd_df <- df %>%
    dplyr::select(map_type, log2FoldChange, taqman_fc) %>%
#    mutate(sq_error = sqrt((log2FoldChange - taqman_fc)^2)) %>%
    mutate(sq_error = log2FoldChange - taqman_fc) %>%
    ungroup()  %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) 

pval <- kruskal.test(sq_error~factor(map_type), data=rmsd_df)$p.value

colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
rmse_bar <- ggplot(data=rmsd_df, aes(x = map_type, y = (sq_error), fill = map_type)) +
    geom_violin() +
    labs(x = 'Pipeline', y = expression(paste(Delta~'log2(fold change)'))) +
    theme(legend.position = 'none') +
    scale_fill_manual(values=colors)

p <- plot_grid(rmse_bar, gene_roc_p, 
               labels = letters[1:2], 
               label_size=20,ncol=1,
               align='v', axis='l')
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/taqman_roc.pdf')
save_plot(p, file=figurename,  base_width=10, base_height=8) 
message('Written: ', figurename)
    
missing_gene <- df %>% 
#    filter(comparison == 'AB') %>% 
    filter(comparison == 'CD') %>% 
    select(-pvalue, - real_FC, - log2FoldChange) %>% 
    spread(map_type, taqman_fc) %>% 
    filter(is.na(`HISAT2+featureCounts`))

df %>% group_by(comparison, map_type) %>% summarize(n())

maqc_roc_table <- str_c(figurepath, '/maqc_roc.csv')
roc_plot_df %>% 
    mutate(map_type = str_replace(map_type,' \\(AUC: 0.[0-9]+\\)','')) %>%
    select(map_type, auc) %>%
    rename(auc_maqc = auc) %>%
    distinct() %>%
    write_csv(maqc_roc_table)
