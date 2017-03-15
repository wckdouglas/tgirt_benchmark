#!/usr/bin/env Rscript

library(robustbase)
library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(broom)
library(purrr)
library(cowplot)
library(feather)

# read gene table
gene_file <- '/stor/work/Lambowitz/ref/GRCh38/transcripts.tsv' %>%
    read_tsv()  %>%
    dplyr::rename(target_id=t_id) %>%
    tbl_df

# read all tables ====================================================
project_path <- '/stor/scratch/Lambowitz/cdw2854/bench_marking'
alignment_free_df <- project_path %>%
    str_c('/alignment_free/countFiles/sleuth_results.feather') %>%
    read_feather() %>%
    mutate(logFC = -b) %>%
    dplyr::rename(id = gene_id) %>%
    select(id,sample_base,sample_test, name,type,qval, mean_obs, logFC) %>%
    gather(variable, value, -id,-name,-type,-sample_base,-sample_test) %>%
    mutate(variable = case_when(
                    .$variable == 'qval' ~ 'adj.P.Val',
                    .$variable == 'mean_obs' ~ 'AveExpr',
                    .$variable == 'logFC' ~ 'logFC'
    ))  %>%
    mutate(variable = str_c(variable,'_', sample_base, sample_test)) %>%
    select(-sample_base,-sample_test) %>%
    spread(variable, value)  %>%
    select(-name,-type) %>%
    mutate(map_type = 'Kallisto')

genome_df <- project_path %>%
    str_c('/genome_mapping/deseq_genome.feather') %>%
    read_feather() %>%
    tbl_df

#z <- 1.43`
fc_df <- genome_df %>%
    rbind(alignment_free_df) %>%
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
    mutate(analytic_type = ifelse(grepl('multimap',map_type),'Genome Mapping','Alignment-free')) %>%
    tbl_df
    
quantil_table <- fc_df %>% 
    group_by(map_type, analytic_type) %>% 
    do(tidy(t(quantile(.$AveExpr_AB, c(.99,.9,.75)))) ) %>%
    ungroup()

fc_df <- fc_df %>% 
    inner_join(quantil_table) %>% 
    mutate(labeling = case_when(
                        .$AveExpr_AB > .$`X99.` ~ 'Top 1%', 
                        .$AveExpr_AB > .$`X90.` ~ 'Top 10%',
                        .$AveExpr_AB > .$`X75.` ~ 'Top 25%',
                        .$AveExpr_AB <= .$`X75.`~ 'Others'))  %>%
    tbl_df

rsquare <- fc_df %>%
    group_by(map_type,analytic_type)%>%#,slope) %>%
    summarize(rs = 1 - sum(error^2)/sum((logFC_CD - mean(logFC_CD))^2),
              samplesize=n()) %>%
    mutate(ypos = 1:(nrow(.)/2))

colors <- c('gray','salmon','light blue', 'goldenrod1')
p <- ggplot() + 
    geom_point(data = fc_df %>% arrange(AveExpr_AB),# %>% filter(grepl('^TR|tR|[ACTG]{3}$',id)), 
               aes(color = labeling,x=logFC_AB, y = logFC_CD, alpha=labeling))  + 
    labs(x = 'log(fold change AB)', y = 'log(fold change CD)', color = ' ') +
    facet_grid(.~analytic_type+map_type) +
    xlim(-10,10) +
    ylim(-2,2) +
    geom_text(x = 7, y = -1.5, data = rsquare, 
              aes(label = paste0('R^2: ',signif(rs,3))), parse=T) +
    geom_text(x = 7, y = -1.75, data = rsquare, 
              aes(label = paste0('n = ',samplesize))) +
    geom_line(data = fc_df, aes(x = logFC_AB, y = predict), color ='red') +
    geom_line(data = fc_df, aes(x = logFC_AB, y = upper_error), color ='red') +
    geom_line(data = fc_df, aes(x = logFC_AB, y = lower_error), color ='red') +
    scale_color_manual(values = colors) +
    scale_alpha_manual(values = c(0.1,1,1,1), guide=F) +
    panel_border()


# make type FC
#======set up change type funciton
ncRNA=c("sense_intronic","3prime_overlapping_ncrna",'processed_transcript',
        'sense_overlapping','Other_lncRNA')
smncRNA=c('misc_RNA','snRNA','piRNA')
large_rRNA=c('28S_rRNA','18S_rRNA')
small_rRNA=c('rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA')
protein_coding = c('protein_coding','TR','IG')

changeType <- function(x){
    if (x %in% ncRNA){
        'Other ncRNA'
    }else if (grepl('TR|IG|protein',x)){
        'Protein coding'
    }else if (grepl('Mt_',x)){
        'Mt'
    }else if (grepl('tRNA',x)){
        'tRNA'
    }else if (x %in% small_rRNA){
        '5/5.8S rRNA'
    }else if (x %in% large_rRNA){
        '18/28S rRNA'
    }else if (x %in% smncRNA){
        'Other sncRNA'
    }else if (grepl('pseudogene',x)){
        'Pseudogenes'
    }else {
        x
    }
}

fc_type_df <- gene_file %>% 
    dplyr::rename(id=gene_id) %>% 
    select(id,name,type) %>% 
    unique %>% 
    inner_join(fc_df) %>% 
    mutate(type = sapply(type,changeType))

colors_type <- RColorBrewer::brewer.pal(12,'Paired')
type_p <- ggplot() + 
    geom_point(data = fc_type_df %>% arrange(AveExpr_AB), 
               aes(shape = labeling, color = type,x=logFC_AB, y =logFC_CD, alpha=labeling))  + 
    labs(x = 'log(fold change AB)', y = 'log(fold change CD)', color = ' ', shape = ' ') +
    facet_grid(.~analytic_type+map_type) +
    xlim(-10,10) +
    ylim(-2,2) +
    geom_text(x = 7, y = -1.5, data = rsquare, 
              aes(label = paste0('R^2: ',signif(rs,3))), parse=T) +
    geom_text(x = 7, y = -1.75, data = rsquare, 
              aes(label = paste0('n = ',samplesize))) +
    geom_line(data = fc_df, aes(x = logFC_AB, y = predict), color ='red') +
    scale_color_manual(values = colors_type) +
    scale_alpha_manual(values = c(0.1,1,1,1), guide=F) +
    panel_border()
p <- plot_grid(p,type_p,ncol=1,align='v')
figurepath <- str_c(project_path, '/figures')
figurename <- str_c(figurepath, '/fold_change.png')
save_plot(p , file=figurename,  base_width=14, base_height=9) 
message('Saved: ', figurename)