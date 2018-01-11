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
    dplyr::select(gene_id, name,type) %>% 
    dplyr::rename(id = gene_id) %>%
    mutate(id = ifelse(type=='tRNA', name, id)) %>%
    mutate(id = ifelse(grepl('MT',id), str_replace(id,'[0-9]$',''), id)) %>%
    distinct() %>%
#    mutate(type= ifelse(grepl('MT-T',name),'tRNA',type)) %>% 
    tbl_df

gene_length <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/genes.length' %>%
    read_tsv() 

# check <- gene_file %>% group_by(id) %>% summarize(a=n()) %>% filter(a>1)
# read all tables ====================================================
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
df <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern = '.feather', full.names=T) %>%
    .[!grepl('abundance|tpm|_[0-9]+|fc_table',.)] %>%
    map_df(read_feather) %>%
    mutate(map_type = case_when(
                grepl('Conventional',.$map_type) ~ "HISAT2+featureCounts",
                grepl('Customized', .$map_type) ~ "TGIRT-map",
                TRUE~ .$map_type)) %>%
    tbl_df
    
    
#z <- 1.43`
fc_df <- df %>%
    mutate(comparison = str_replace(comparison, ' vs ','')) %>%
    gather(variable, value, -id, -comparison,-map_type) %>%
    filter(!grepl('lfcSE|stat', variable)) %>%
    mutate(variable = str_c(variable, comparison, sep='_')) %>%
    dplyr::select(-comparison) %>%
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
    dplyr::select(-z,-k1,-k2) %>%
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
                   dplyr::select(id,name,type) %>% 
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
                        .$baseMean_AB <= .$`X75.`~ 'Bottom 75%'))  %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) %>%
    tbl_df

rsquare <- fc_df %>%
    ungroup()%>%
    group_by(map_type,analytic_type, labeling)%>%#,slope) %>%
    summarize(
        rs = caret::R2(predict, logFC_CD, formula='traditional'),
        sum_err = sum(error^2),
        sum_var = sum((logFC_CD - mean(logFC_CD))^2),
        samplesize=n()
    ) %>%
    mutate(ypos = 1) %>%
    mutate(ypos = cumsum(ypos)) %>%
    ungroup() %>%
    mutate(rs = formatC(round(rs,3), format='f',digits=3)) %>%
    tbl_df



fc_type_df <- gene_file %>% 
    dplyr::select(id,name,type) %>% 
    distinct() %>% 
    inner_join(fc_df)  %>%
    tbl_df

rmse_fc <- fc_type_df %>% 
    group_by(map_type) %>% 
    summarize(
        rs = caret::R2(predict, logFC_CD, formula='traditional'),
        sum_err = sum(error^2),
        sum_var = sum((logFC_CD - mean(logFC_CD))^2),
        samplesize=n()
    ) %>%
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
    dplyr::select(map_type, padj_AB, id) %>%
    unique() %>%
    spread(map_type, padj_AB) %>%
    tbl_df
sig_scatter <- GGally::ggpairs(sig_DE_df, 2:5)
dev.off()
message('Saved: ', figurename)


per_type_r2_df <- fc_type_df  %>% 
    group_by(type,map_type,analytic_type) %>%     
    summarize(
            rmse = caret::RMSE(predict, logFC_CD),
            rs = caret::R2(predict, logFC_CD, formula='traditional'),
            rs_Cor = caret::R2(predict, logFC_CD),
            samplesize=n()) %>%
    ungroup() %>%
    mutate(prep = str_c(map_type,sep='\n')) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',prep),1,2)) %>%
    mutate(prep = forcats::fct_reorder(prep, pipeline_type)) %>%
    tbl_df

colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
error_table <- str_c(figurepath, '/pertype_rmse.csv')
per_type_r2_df %>% 
    mutate(value = str_c('\\shortstack{',signif(rmse,3),'\\\\ (',samplesize,')}')) %>%
    dplyr::select(type, map_type, value)  %>%
    spread(map_type, value) %>%
    mutate(type = str_c('\\textbf{',type,'}')) %>%
    dplyr::rename(`RNA type` = type) %>%
    write_csv(error_table)

rmse_type_p <- ggplot(data=per_type_r2_df %>% filter(!grepl('rRNA',type)), 
       aes(x=prep,y=rmse, fill=prep)) + 
    geom_bar(stat='identity') +
    #geom_text(aes(label = str_c(samplesize,'\n(',signif(rmse,3),')')), size=4, hjust=0,vjust=0.5,angle= 90)+
    facet_wrap(~type) +
    labs(color = ' ', fill= ' ', x= 'Pipeline', y = 'Root mean square error', parse=T) +
    panel_border()+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position = c(0.65,0.2)) +
    ylim(0,0.8) + 
    scale_fill_manual(values=colors)

r2_type_p <- ggplot(data=per_type_r2_df %>% filter(!grepl('rRNA',type)), 
                      aes(x=prep,y=rs, fill=prep)) + 
    geom_bar(stat='identity') +
    #geom_text(aes(label = str_c(samplesize,'\n(',signif(rmse,3),')')), size=4, hjust=0,vjust=0.5,angle= 90)+
    facet_wrap(~type, scale='free') +
    labs(color = ' ', fill= ' ', x= 'Pipeline', y = expression('R^2'), parse=T) +
    panel_border()+
    theme(axis.text.x=element_blank())+
    theme(axis.ticks.x=element_blank())+
    theme(legend.key.height = unit(2,'line')) +
    theme(legend.position = c(0.65,0.2)) +
#    ylim(0,0.8) + 
    scale_fill_manual(values=colors) +
    geom_hline(yintercept = 0)


colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
tRNA_error <- fc_type_df %>% 
    filter(type=='tRNA') %>%
    dplyr::select(map_type,error, name) %>%
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
        theme(legend.position = 'none')+
        scale_color_manual(values=colors) +
        draw_text('Alignment-free', 24, 5) +
        draw_text('Alignment-based', 53, 2)
    
font_size = 12
p <- plot_grid(rmse_type_p+
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size),
                         strip.text = element_text(size=font_size)), 
               tRNA_error + 
                   theme(axis.text = element_text(size=font_size),
                         axis.title = element_text(size=font_size)),
               ncol=1,align='v', axis='l',
               labels = letters[1:2], label_size=20,
               rel_heights=c(1.5,1))
figurename <- str_c(figurepath, '/fold_change_all_gene.pdf')
save_plot(p , file=figurename,  base_width=9, base_height=11) 
message('Saved: ', figurename)


#p2 <- plot_grid(type_p + 
#                   theme(axis.text = element_text(size=font_size),
#                         axis.title = element_text(size=font_size),
#                         strip.text = element_text(size=font_size)),
#                s_p +
#                    ggrepel::geom_label_repel(data = fc_df %>% 
#                                  filter(labeling %in% c('Top 1%','Top 10%')) %>% 
#                                  group_by(labeling, analytic_type, map_type) %>%
#                                  top_n(5,error),
#                                  alpha=1,
#                              aes(x=logFC_AB, y = logFC_CD, label = name, color = labeling),
#                              show.legend = F),
#               ncol=1, align='v', axis='l', 
#               labels=letters[1:2], label_size=30)

colors <- c("#E69F00", "#F0E442", "#009E73", "#56B4E9")
p2 <- rsquare %>% 
    dplyr::select(map_type, labeling, rs) %>% 
    rbind(rmse_fc %>% 
              dplyr::select(map_type, rs) %>% 
              mutate(labeling = 'Total RNA')
          ) %>%
    mutate(labeling = factor(labeling, levels=c('Total RNA','Top 1%', 'Top 10%', 'Top 25%','Bottom 75%'))) %>%
    mutate(rs = as.numeric(rs)) %>%
    ggplot(aes(x=map_type, y=rs, fill=map_type)) +
        geom_bar(stat='identity') + 
        facet_grid(~labeling) +
        scale_fill_manual(values=colors) +
        labs(x='Pipeline', y = expression(R^2), fill= ' ') +
        theme(axis.text.x=element_blank()) +
        ylim(0,1)



fc_type_df %>% 
    filter(type=='tRNA') %>% 
    dplyr::select(id,name, type, map_type, baseMean_AB) %>% 
    spread(map_type, baseMean_AB) %>% 
    filter(is.na(Salmon))

fc_type_df %>% 
    filter(type=='vaultRNA') %>% 
    dplyr::select(id,name, type, map_type, baseMean_AB) %>% 
    spread(map_type, baseMean_AB) 

rmse_table <- str_c(figurepath, '/rmse_fc.csv')
per_type_r2_df %>% 
    dplyr::select(rmse,map_type,type,analytic_type) %>%
    spread(type, rmse) %>%
    write_csv(rmse_table)

fc_type_df %>% 
        filter(type == 'Other sncRNA') %>% 
        ggplot(aes(x = predict, y = logFC_CD)) + 
            stat_density2d(geom='raster', contour=F, aes(fill=..density..)) + 
            viridis::scale_fill_viridis() + 
            facet_grid(.~map_type) + 
            geom_abline(color='red', slope=1, intercept= 0) + 
            geom_smooth(method='lm', color= 'green', alpha = 0.6)


quantile_length_r2 <- fc_df %>% 
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
        rs_cor = caret::R2(predict, logFC_CD),
        rmse = caret::RMSE(predict, logFC_CD),
        rmse = sqrt(mean(error ^ 2)),
        uFC = mean(logFC_CD),
        varFC = var(logFC_CD),
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
    
length_tile_p <- ggplot(data=quantile_length_r2, 
                        aes(x=map_type, y = rs, fill = map_type)) +
    geom_bar(stat='identity') +
    labs(x = 'Pipeline', y = expression('R^2'))+
    facet_grid(~length_tile_group) +
    scale_fill_manual(values=colors) +
    labs(x='Pipeline', y = expression(R^2), fill= ' ') +
    theme(axis.text.x=element_blank()) 

r2_p <- plot_grid(p2, length_tile_p, ncol=1, align = 'v',
                  labels = c('a','b'), axis='l')
figurename <- str_c(figurepath, '/fold_change_supplementary.png')
save_plot(r2_p , file=figurename,  base_width=10, base_height=10) 
message('Saved: ', figurename)


quantile_length <- fc_df %>% 
    inner_join(gene_length)  %>%
#    write_feather('/stor/work/Lambowitz/cdw2854/bench_marking/DEgenes/fc_table.feather')
    group_by(map_type) %>%
    nest() %>%
    mutate(data = map(data, ~mutate(., length_tile = ntile(gene_length, 4)))) %>%
    unnest() %>%
    mutate(length_tile_group = case_when(
                        .$length_tile == 4 ~ '0-25%\n(Longest genes)',
                        .$length_tile == 3 ~ '25-50%',
                        .$length_tile == 2 ~ '50-75%',
                        .$length_tile == 1 ~ '75-100%\n(Shortest genes)'
    )) %>% 
    filter(map_type == 'Kallisto') %>%
    tbl_df
    

r2_df <- quantile_length_r2 %>% 
    mutate(r2_label = round(rs,3)) %>%
    mutate(r2_label = paste('R^2','==',r2_label)) %>%
    mutate(rmse_label = round(rmse,3)) %>%
    mutate(rmse_label = paste('RMSE','==',rmse_label)) %>%
    filter(map_type == 'Kallisto') %>%
    dplyr::select(r2_label, rmse_label, length_tile_group, map_type, rs)
colors <- RColorBrewer::brewer.pal(8, 'Dark2')
s_p <- ggplot() + 
    stat_density2d(geom = "raster", 
                    data = quantile_length %>% arrange(-length_tile),
                    aes(x=predict, y = logFC_CD, fill = ..density..), 
                    contour = FALSE) +
#    geom_smooth(data = quantile_length %>% arrange(-length_tile),
#                aes(x=predict, y = logFC_CD), method='lm', 
#                color='purple', alpha=0.5)+
    facet_wrap(~length_tile_group, scale='free') +
    labs(x = 'Expected fold change CD (log2)', y = 'Measured fold change CD (log2)', 
         color = ' ', fill = 'Density') +
    xlim(-1,1) +
    ylim(-1,1) +
    geom_abline(intercept = 0, slope = 1, color = 'red', alpha=0.6) +
    geom_text(data=r2_df, x = -0.4, y = 0.8,
              aes(label = r2_label), face='bold', size = 7, 
              color='white', parse=T)+
    geom_text(data=r2_df, x = -0.4, y = 0.6,
              aes(label = rmse_label), face='bold', size = 7, 
              color='white', parse=T)+
    viridis::scale_fill_viridis() +
    panel_border()  +
    theme(strip.text = element_text(size = 16))
message('Plotted S curve')
figurename <- str_c(figurepath, '/kallisto_r2_splot.png')
save_plot(s_p, file=figurename,  base_width=8, base_height=8) 
message('Saved: ', figurename)

# 
# abund_r2_df <- rsquare %>% 
#     mutate(r2_label = paste('R^2','==',rs)) %>%
#     dplyr::select(r2_label, labeling, map_type, rs) %>%
#     mutate(labeling = factor(labeling, 
#                              levels = rev(c('Bottom 75%','Top 25%','Top 10%','Top 1%'))))
# abund_s_p <- ggplot(data= fc_df %>% 
#                         mutate(labeling = factor(labeling, 
#                             levels = rev(c('Bottom 75%','Top 25%','Top 10%','Top 1%'))))) + 
#     stat_density2d(geom = "raster", 
#                    #aes(x=logFC_AB, y = logFC_CD, fill = ..density..), 
#                    aes(x=predict, y = logFC_CD, fill = ..density..), 
#                    contour = FALSE) +
#     facet_grid(labeling~map_type, scale='free') +
#     labs(x = 'Expected fold change AB (log2)', y = 'Measured fold change CD (log2)', 
#          color = ' ', fill = 'Density') +
#     #    facet_grid(.~analytic_type+map_type) +
#     xlim(-2.5,2.5) +
#     ylim(-2,2) +
#     #    geom_text(x = -7, data = rsquare, alpha=1,
#     #              aes(y=2-ypos/5, 
#     #                  label = paste0('R^2: ',as.character(rs)), 
#     #                  color = labeling), parse=T, show.legend=F) +
#     #    geom_text(x = 7, data = rsquare, alpha=1,
#     #              aes(y=-1-ypos/5,label = paste0('n = ',samplesize), color = labeling), show.legend = F) +
#     #geom_line(aes(x = logFC_AB, y = predict), 
#     #          color ='red', alpha=0.6) +
#     geom_abline(intercept = 0, slope = 1, color = 'red', alpha=0.6) +
#     geom_text(data=abund_r2_df, x = -1, y = 1.5,
#               aes(label = r2_label),
#               color='white', parse=T)+
#     viridis::scale_fill_viridis() +
#     panel_border() 
# 
# supp_s_p <- plot_grid(abund_s_p, s_p, 
#                       labels = c('a','b'), 
#                       ncol=1)
# figurename <- str_c(figurepath, '/fold_change_r2_splot.png')
# save_plot(supp_s_p, file=figurename,  base_width=10, base_height=13) 
# message('Saved: ', figurename)
