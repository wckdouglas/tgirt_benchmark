#!/usr/bin/env Rscript

#library(robustbase)
library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(broom)
library(purrr)
library(tidyr)
library(cowplot)
library(feather)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv'  %>%
    read_tsv() %>%
    select(gene_id, name,type) %>% 
    dplyr::rename(id = gene_id) %>%
    mutate(id = str_replace(id, '_gene$','')) %>%
    mutate(id = ifelse(type=='tRNA', name, id)) %>%
    mutate(id = ifelse(grepl('MT',id), str_replace(id,'[0-9]$',''), id)) %>%
    distinct() %>%
    distinct() %>%
#    mutate(type= ifelse(grepl('MT-T',name),'tRNA',type)) %>% 
    tbl_df

# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|Conven|kalli|fc_table',.)] %>%
    map_df(read_feather) %>%
    mutate(map_type = case_when(
                        grepl('align', .$map_type)~'Salmon: Alignment-based',
                        grepl('_[0-9]+', .$map_type)~str_replace(.$map_type,'_','\nkmer='),
                        grepl('Salmon', .$map_type)~'Salmon\nkmer=31',
                        grepl('Cust',.$map_type) ~ "TGIRT-map",
                        TRUE~.$map_type))
    
    
#z <- 1.43`
fc_df <- df %>%
    mutate(comparison = str_replace(comparison, ' vs ','')) %>%
    gather(variable, value, -id, -comparison,-map_type) %>%
    filter(!grepl('lfcSE|stat', variable)) %>%
    mutate(variable = str_c(variable, comparison, sep='_')) %>%
    select(-comparison) %>%
    spread(variable, value) %>%
    dplyr::rename(logFC_AB = log2FoldChange_AB) %>%
    dplyr::rename(logFC_CD = log2FoldChange_CD) %>%
    filter(!is.na(logFC_AB)) %>%
    filter(!is.na(logFC_CD)) %>% 
    mutate(AB = 2^logFC_AB) %>%
    mutate(z = 1.43) %>%
    mutate(k1 = 3*z/(3*z+1)) %>%
    mutate(k2 = z/(z+3)) %>%
    mutate(predict = log2(k1 * 2^(logFC_AB) + (1-k1)) - log2(k2 * 2^(logFC_AB) + (1-k2))) %>%
    mutate(z = 1.43*1.1) %>%
    mutate(k1 = 3*z/(3*z+1)) %>%
    mutate(k2 = z/(z+3)) %>%
    mutate(upper_error = log2(k1 * 2^(logFC_AB) + (1-k1)) - log2(k2 * 2^(logFC_AB) + (1-k2))) %>%
    mutate(z = 1.43*0.9) %>%
    mutate(k1 = 3*z/(3*z+1)) %>%
    mutate(k2 = z/(z+3)) %>%
    mutate(lower_error = log2(k1 * 2^(logFC_AB) + (1-k1)) - log2(k2 * 2^(logFC_AB) + (1-k2))) %>%
    select(-z,-k1,-k2) %>%
    mutate(error = logFC_CD - predict) %>% 
    ungroup() %>%
    mutate(analytic_type = ifelse(grepl('pipe',map_type),'Genome Mapping','Alignment-free')) %>%
    inner_join(gene_file %>% 
                   select(id,name,type) %>% 
                   unique,
               by = 'id') %>%
    mutate(map_type = str_replace(map_type, 'kmer','k-mer')) %>%
    tbl_df
    
quantil_table <- fc_df %>% 
    group_by(map_type, analytic_type) %>% 
    do(tidy(t(quantile(.$baseMean_AB, c(.99,.9,.75)))) ) %>%
    ungroup()

fc_df <- fc_df %>% 
    inner_join(quantil_table, by = c("map_type", "analytic_type")) %>% 
    mutate(labeling = case_when(
                        .$baseMean_AB > .$`X99.` ~ 'Top 1%', 
                        .$baseMean_AB > .$`X90.` ~ 'Top 10%',
                        .$baseMean_AB > .$`X75.` ~ 'Top 25%',
                        .$baseMean_AB <= .$`X75.`~ 'Others'))  %>%
    tbl_df

rsquare <- fc_df %>%
    ungroup()%>%
    group_by(map_type,analytic_type, labeling)%>%#,slope) %>%
    summarize(
        sum_err = sum(error^2),
        sum_var = sum((logFC_CD - mean(logFC_CD))^2),
        samplesize=n()
    ) %>%
    mutate(rs =  1 - sum_err/sum_var) %>%
    mutate(ypos = 1) %>%
    mutate(ypos = cumsum(ypos)) %>%
    ungroup() %>%
    mutate(rs = formatC(round(rs,3), format='f',digits=3)) %>%
    tbl_df

fc_type_df <- gene_file %>% 
    select(id,name,type) %>% 
    distinct() %>% 
    inner_join(fc_df) 

rmse_fc <- fc_type_df %>% 
    group_by(map_type) %>% 
    summarize(
        sum_err = sum(error^2),
        sum_var = sum((logFC_CD - mean(logFC_CD))^2),
        samplesize=n()
    ) %>%
    mutate(rs =  1 - sum_err/sum_var) %>%
    mutate(rs = signif(rs, 3)) %>%
    tbl_df


per_type_r2_df <- fc_type_df  %>% 
    group_by(type,map_type,analytic_type) %>%     
    summarize(
            rmse = sqrt(mean(error^2)),
            samplesize=n()) %>%
    ungroup() %>%
    mutate(prep = str_c(map_type,sep='\n')) %>%
    tbl_df

error_table <- str_c(figurepath, '/pertype_rmse_kmer.csv')
per_type_r2_df %>% 
#    mutate(value = str_c('\\shortstack{',signif(rmse,3),'\\\\ (',samplesize,')}')) %>%
    mutate(value = str_c('= "',signif(rmse,3),'" & CHAR(13) & "(',samplesize,')"')) %>%
    select(type, map_type, value)  %>%
    spread(map_type, value) %>%
#    mutate(type = str_c('\\textbf{',type,'}')) %>%
    dplyr::rename(`RNA type` = type) %>%
    select(-`TGIRT-map`) %>%
    write_csv(error_table)


colors <- c('#BE8088', '#B88573', '#9C915A', '#789A63','#699695', '#F0E442')
rmse_type_p <- ggplot(data=per_type_r2_df %>% filter(!grepl('rRNA',type)) %>% mutate(prep = str_replace(prep,'\n',' ')), 
       aes(x=prep,y=rmse, fill=prep)) + 
    geom_bar(stat='identity') +
#    geom_text(aes(label = str_c(samplesize,'\n(',signif(rmse,3),')')), size=4, hjust=0,vjust=0.5,angle= 90)+
    facet_wrap(~type) +
    labs(color = ' ', fill= ' ', x= 'Pipeline', y = 'Root mean square error', parse=T) +
    panel_border()+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(legend.key.height = unit(1,'line')) +
    theme(legend.position = c(0.6,0.15)) +
    scale_fill_manual(values=colors) 

tRNA_error <- fc_type_df %>% 
    filter(type=='tRNA') %>%
    select(map_type,error, name) %>%
    mutate(error = abs(error)) %>%
    group_by(map_type) %>%
    nest() %>%
    mutate(data = map(data, function(d) d %>% 
                          arrange(error) %>% 
                          mutate(cum_error = cumsum(error)) %>%
                          mutate(id = 1:nrow(.)))) %>%
    unnest() %>%
    ggplot(aes(x=id, y=cum_error)) +
        geom_line(size=2, aes(color=map_type)) +
        labs(y = 'Cumulative error\n(log2 fold change)',
             x='Number of tRNA', color = ' ') +
        scale_x_continuous(breaks=seq(0,56,8), limits=c(1,56)) +
        scale_color_manual(values=colors) +
        theme(legend.position = 'none') +
        draw_text('k = 11 or 31', 40, 13) +
        draw_text('k = 15 or 21', 53, 11) +
        draw_text('TGIRT-map', 50, 2) 


quantile_length <- fc_df %>% 
    #    filter(type != 'Other ncRNA') %>%
    inner_join(gene_length)  %>%
    group_by(map_type) %>%
    nest() %>%
    mutate(data = map(data, ~mutate(., length_tile = ntile(gene_length, 4)))) %>%
    unnest(data) %>%
    ungroup() %>%
    group_by(map_type, length_tile) %>%
    summarize(
        rs = caret::R2(predict, logFC_CD, formula='traditional'),
        rmse = sqrt(mean(error ^ 2)),
        samplesize=n()
    ) %>%
    ungroup %>%
    mutate(length_tile_group = case_when(
        .$length_tile == 4 ~ '0-25%\n(Longest genes)',
        .$length_tile == 3 ~ '25-50%',
        .$length_tile == 2 ~ '50-75%',
        .$length_tile == 1 ~ '75-100%\n(Shortest genes)'
    )) %>%
    tbl_df

length_tile_p <- ggplot(data=quantile_length, aes(x=map_type, y = rs, fill = map_type)) +
    geom_bar(stat='identity') +
    labs(x = 'Pipeline', y = expression('R^2'))+
    facet_grid(~length_tile_group) +
    scale_fill_manual(values=colors) +
    labs(x='Pipeline', y = expression(R^2), fill= ' ') +
    theme(axis.text.x=element_blank())  +
    theme(legend.key.height = unit(2,'line')) 

p <- plot_grid(length_tile_p+
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size),
                         strip.text = element_text(size=font_size)), 
               tRNA_error + 
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size)),
               ncol=1,align='v', axis='l',
               labels = letters[1:2], label_size=20,
               rel_heights=c(1.5,1))
figurename <- str_c(figurepath, '/fold_change_all_gene_kmer.png')
save_plot(p, file=figurename,  base_width=8, base_height=8) 
message('Saved: ', figurename)