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
gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/all_genes.tsv'  %>%
    read_tsv() %>%
    select(gene_id, name,type) %>% 
    dplyr::rename(id = gene_id) %>%
    mutate(id = ifelse(type=='tRNA', name, id)) %>%
    mutate(id = ifelse(grepl('MT',id), str_replace(id,'[0-9]$',''), id)) %>%
    distinct() %>%
#    mutate(type= ifelse(grepl('MT-T',name),'tRNA',type)) %>% 
    tbl_df

# check <- gene_file %>% group_by(id) %>% summarize(a=n()) %>% filter(a>1)
# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|_[0-9]+',.)] %>%
    map_df(read_feather) %>%
    mutate(map_type = case_when(
                grepl('Conventional',.$map_type) ~ "HISAT2+FeatureCounts",
                grepl('Customized', .$map_type) ~ "TGIRT-map",
                TRUE~ .$map_type)) %>%
    tbl_df

 
    
    
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
    #unnest(data) %>%
    #group_by(map_type) %>%
    #nest() %>%
    #mutate(model = map(data, ~ nlrob(logFC_CD ~ log2((3*z/(3*z+1)) * AB + (1- (3*z/(3*z+1)))) - log2((z/(3+z)) * AB+(1-(z/(3+z)))),
    #                                 data=., start=list(z=1.43)))) %>%
    #mutate(predict = map2(model, data, predict)) %>%
    #unnest(data,predict) %>%
    mutate(error = logFC_CD - predict) %>% 
    ungroup() %>%
    mutate(analytic_type = ifelse(grepl('HI|TG',map_type),'Genome Mapping','Alignment-free')) %>%
    inner_join(gene_file %>% 
                   select(id,name,type) %>% 
                   unique,
               by = 'id') %>%
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

colors <- c('gray','salmon','light blue', 'goldenrod1')
s_p <- ggplot() + 
    geom_point(data = fc_df %>% arrange(baseMean_AB),# %>% filter(grepl('^TR|tR|[ACTG]{3}$',id)), 
               aes(color = labeling,x=logFC_AB, y = logFC_CD, alpha=labeling))  + 
    labs(x = 'Fold change AB (log2)', y = 'Fold change CD (log2)', color = ' ') +
#    facet_grid(.~analytic_type+map_type) +
    facet_grid(.~map_type) +
    xlim(-10,10) +
    ylim(-2,2) +
    geom_text(x = -7, data = rsquare, alpha=1,
              aes(y=2-ypos/5, 
                  label = paste0('R^2: ',as.character(rs)), 
                  color = labeling), parse=T) +
    geom_text(x = 7, data = rsquare, alpha=1,
              aes(y=0-ypos/5,label = paste0('n = ',samplesize), color = labeling)) +
    geom_line(data = fc_df, aes(x = logFC_AB, y = predict), color ='red') +
    geom_line(data = fc_df, aes(x = logFC_AB, y = upper_error), color ='red') +
    geom_line(data = fc_df, aes(x = logFC_AB, y = lower_error), color ='red') +
    scale_color_manual(values = colors) +
    scale_alpha_manual(values = c(0.1,1,1,1), guide=F) +
    panel_border()
message('Plotted S curve')


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


colors_type <- RColorBrewer::brewer.pal(12,'Paired')
colors_type <- c(colors_type, 'darkgrey')
type_p <- ggplot() + 
    geom_point(data = fc_type_df %>% 
                    arrange(baseMean_AB),
               aes(shape = labeling, color = type, x=logFC_AB, y =logFC_CD, alpha=labeling))  + 
    geom_text(x = -7, y =2, data = rmse_fc, 
              aes(label = paste0('R^2: ',as.character(rs))), parse=T) +
    geom_text(x = -7, y =1.7, data = rmse_fc, 
              aes(label = paste0('(n=',as.character(samplesize),')'))) +
    labs(x = 'Fold change AB (log2)', y = 'Fold change CD (log2)', color = ' ', shape = ' ') +
#    facet_grid(.~analytic_type+map_type) +
    facet_grid(.~map_type) +
    xlim(-10,10) +
    ylim(-2,2) +
    geom_line(data = fc_df, aes(x = logFC_AB, y = predict), color ='red') +
    scale_color_manual(values = colors_type) +
    scale_alpha_manual(values = c(0.1,1,1,1), guide=F) +
    panel_border()


figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/signif_genes_corr.png')
png(figurename)
sig_DE_df <- fc_df %>% 
    filter(padj_AB<0.05) %>% 
    mutate(map_type =  make.names(str_c(analytic_type, '_',map_type))) %>%
    select(map_type, padj_AB, id) %>%
    unique() %>%
    spread(map_type, padj_AB) %>%
    tbl_df
sig_scatter <- GGally::ggpairs(sig_DE_df, 2:5)
dev.off()
message('Saved: ', figurename)


per_type_r2_df <- fc_type_df  %>% 
    group_by(type,map_type,analytic_type) %>%     
    summarize(
            rmse = sqrt(mean(error^2)),
            samplesize=n()) %>%
    ungroup() %>%
    mutate(prep = str_c(map_type,sep='\n')) %>%
    tbl_df
rmse_type_p <- ggplot(data=per_type_r2_df %>% filter(!grepl('rRNA',type)), 
       aes(x=prep,y=rmse, fill=prep)) + 
    geom_bar(stat='identity') +
    geom_text(aes(label = str_c(samplesize,'\n(',signif(rmse,3),')')), size=5, hjust=0,vjust=0.5,angle= 90)+
    facet_wrap(~type) +
    labs(color = ' ', fill= ' ', x= ' ', y = 'Root mean square error', parse=T) +
    panel_border()+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position = c(0.9,0.2)) +
    ylim(0,0.8)


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
    ggplot(aes(x=id, y=cum_error, color=map_type)) +
        geom_line(size=2) +
        labs(y = 'Cumulative error\n(log2 fold change)',
             x='Number of tRNA', color = ' ') +
        scale_x_continuous(breaks=seq(0,56,8), limits=c(1,56)) +
        theme(legend.position = 'none')
    
    

font_size = 20
p <- plot_grid(type_p + 
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size),
                         strip.text = element_text(size=font_size)),
               rmse_type_p+
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size),
                         strip.text = element_text(size=font_size)), 
               tRNA_error + 
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size)),
               ncol=1,align='v',
               labels = letters[1:3], label_size=30)
figurename <- str_c(figurepath, '/fold_change_all_gene.png')
save_plot(p , file=figurename,  base_width=15, base_height=22) 
message('Saved: ', figurename)


p2 <- plot_grid(s_p +
                    ggrepel::geom_label_repel(data = fc_df %>% 
                                  filter(labeling %in% c('Top 1%','Top 10%')) %>% 
                                  group_by(labeling, analytic_type, map_type) %>%
                                  top_n(5,error),
                                  alpha=1,
                              aes(x=logFC_AB, y = logFC_CD, label = name, color = labeling)),
               ncol=1,align='v')

figurename <- str_c(figurepath, '/fold_change_supplementary.png')
save_plot(p2 , file=figurename,  base_width=14, base_height=7) 
message('Saved: ', figurename)


fc_type_df %>% 
    filter(type=='tRNA') %>% 
    select(id,name, type, map_type, baseMean_AB) %>% 
    spread(map_type, baseMean_AB) %>% 
    filter(is.na(Salmon))

rmse_table <- str_c(figurepath, '/rmse_fc.csv')
per_type_r2_df %>% 
    select(rmse,map_type,type,analytic_type) %>%
    spread(type, rmse) %>%
    write_csv(rmse_table)
