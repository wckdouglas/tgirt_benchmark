#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(stringi)
library(tidyr)
library(feather)
library(purrr)
library(VennDiagram)
library(GGally)
library(UpSetR)

setwd('~/tgirt_benchmark/find_DE/')
source('ercc_correlation.R')
counts_to_tpm <- function(counts, lengths) {
    #https://www.biostars.org/p/171766/
    #
    rate <- counts / lengths
    out_tpm <- rate / sum(rate,na.rm=T) * 1e6
    return(out_tpm)
}

#read files
get_mix <- function(samplename){
    ifelse(grepl('^ref',samplename), str_sub(samplename,4,4), str_sub(samplename,8,8))
}

get_sample_number <- function(samplename){
    id <- ifelse(grepl('^ref',samplename), str_sub(samplename,5,5), str_sub(samplename,10,10))
    return(as.numeric(id))
}

gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/transcripts.tsv' %>%
    read_tsv() %>%
    select(gene_id, name, type) %>%
    dplyr::rename(id = gene_id) %>%
    unique() %>%
    mutate(id = str_replace(id, '_gene$','')) 
    

gene_length_df <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/genes.length' %>%
    read_tsv()  %>%
    mutate(id = str_replace(id, '_gene$','')) %>%
    mutate(id = str_replace(id, '[0-9]+-[0-9]+$','')) %>%
    group_by(id) %>%
    summarize(gene_length = min(gene_length)) %>%
    ungroup() %>%
    tbl_df
 
#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking_new/bench_marking/'
alignment_free <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
    .[!grepl('genome',.)] %>%
    .[!grepl('_[0-9]+|aligned|bias',.)] %>%
    map_df(read_feather)  %>%
    gather(samplename, abundance, -id, -map_type) %>%
    mutate(id = str_replace(id,'_gene$',''))%>%
    tbl_df

# read genome mapping count files
files <- c(file.path(project_path, 'genome_mapping/conventional/counts/feature_counts.tsv'),
           file.path(project_path,'genome_mapping/tgirt_map/Counts/RAW/combined_gene_count.tsv'))
labels <- c('conventional','customized')
genome_df <- map2(files, labels, function(x,y) read_tsv(x) %>% 
                      mutate(map_type=y) %>%
                      set_names(str_replace_all(names(.),'-','_'))) %>%
    purrr::reduce(rbind) %>%
    mutate(id = str_replace(id,'_gene$',''))%>%
    mutate(id = str_replace(id,'[0-9]+-[0-9]+$',''))%>%
    inner_join(gene_length_df) %>%
    gather(samplename, abundance, -id,-map_type, -gene_length) %>%
    mutate(id = ifelse(!grepl('ERCC',id),str_replace(id, '\\-[0-9]+$',''),id)) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]+$','') ,id)) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]-$','') ,id)) %>%
    mutate(id = str_replace_all(id, '\\([-+]\\)','')) %>%
    mutate(id = ifelse(id %in% c('MT-TL','MT-TS'), str_c(id,'1'), id)) %>%
    group_by(id, samplename, map_type) %>%
    summarize(
        abundance = sum(abundance),
        gene_length = mean(gene_length)
    ) %>%
    ungroup() %>%
    group_by(samplename, map_type) %>%
    do(
       data_frame(
            abundance = counts_to_tpm(.$abundance, .$gene_length),
            id = .$id
        )
    ) %>%
    ungroup() %>%
    mutate(map_type = case_when(
                                grepl('conventional',.$map_type) ~ "HISAT2+featureCounts",
                                grepl('customized', .$map_type) ~ "TGIRT-map",
                                TRUE~ .$map_type)
    ) %>%
    tbl_df

genome_df %>% 
    spread(samplename,abundance) %>%
    write_feather(file.path(project_path, 'DEgenes/genome_abundance_tpm.feather'))

merge_df <- rbind(genome_df, alignment_free) %>%
    mutate(id = str_replace(id, '[0-9]-$','')) %>%
    inner_join(gene_file) %>%
    mutate(id = ifelse(type == 'Mt', name, id)) %>%
    group_by(map_type,samplename, id, name,type) %>%
    summarize(abundance = sum(abundance)) %>%
    tbl_df


plot_df <- merge_df %>% 
    filter(abundance > 0.1) %>%
    group_by(type,samplename, map_type) %>%
    summarize(gene_count = n()) %>%
    ungroup() %>%
    mutate(sample_id = str_extract(samplename, '_([ABCD])_')) %>%
    mutate(sample_id = str_replace_all(sample_id, '_','')) %>%
    mutate(samplenum = str_extract(samplename, '_([1-3])_')) %>%
    mutate(samplenum = str_replace_all(samplenum, '_','')) %>%
    mutate(samplename = str_c(sample_id, samplenum)) %>%
    #mutate(samplename = sample_id) %>%
    dplyr::select(-sample_id, -samplenum ) %>%
    mutate(type = forcats::fct_reorder(type, gene_count , sum)) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type))

friedman_p <-  plot_df %>% 
    group_by(map_type, samplename) %>% 
    summarize(gc = sum(gene_count)) %>% 
    friedman.test(gc~map_type|samplename, data=.) %>%
    .$p.value

#friedman_per_type_p <- plot_df %>%
#    group_by(type) %>%
#    nest() %>%
 #   mutate(friedman = map(data, ~friedman.test(gene_count~map_type|samplename, data=.))) %>%
 #   mutate(p = map_dbl(friedman, function(x) x$p.value)) %>%
 #   unnest(p)

#per type p value
all_comparison <- gtools::permutations(n=4,r=2,
                                      v=unique(as.character(plot_df$map_type)),
                                      repeats.allowed=F)
get_all_p <-function(x1, x2, all_comparison, plot_df){
    compare <- str_c(x1, x2, sep='|')

    wilcox_p <- plot_df %>%
        filter(map_type == x1 | map_type == x2) %>%
        group_by(map_type, samplename) %>%
        summarize(gc = sum(gene_count)) %>% 
        tbl_df
    
    a <- wilcox_p %>% filter(map_type == x1)
    b <- wilcox_p %>% filter(map_type == x2)
    cohen <- effsize::cohen.d(a$gc, b$gc, paired=T) %>%
        .$estimate
    return(data.frame(cohen, compare))
}
pval_df <- all_comparison %>%
    data.frame() %>%
    mutate(wilcox_p = map2(X1, X2, get_all_p, all_comparison, plot_df)) %>%
    unnest(wilcox_p) %>%
    tbl_df

#per type p value
all_comparison = gtools::permutations(n=4,r=2,
                                      v=unique(as.character(plot_df$map_type)),
                                      repeats.allowed=F)
get_p <-function(i, all_comparison, plot_df){
    compare <- str_c(all_comparison[i,1],all_comparison[i,2],sep='|')
    message('Comparing: ', all_comparison[i,1], ' and ', all_comparison[i,2])

    wilcox_per_type_p <- plot_df %>%
        filter(map_type == all_comparison[i,1] | map_type == all_comparison[i,2]) %>%
        mutate(map_type = factor(map_type)) %>%
        arrange(samplename) %>%
        group_by(type) %>%
        nest() %>%
        mutate(wilcox = map(data, ~wilcox.test(gene_count~map_type, data=., paired=T))) %>%
        mutate(p = map_dbl(wilcox, function(x) x$p.value)) %>%
        unnest(p) %>%
        mutate(comparison = compare) %>%
        select(-wilcox, - data)
    return(wilcox_per_type_p)
}
pval_df <- map_df(1:nrow(all_comparison), get_p, all_comparison, plot_df) %>%
    mutate(comparison = str_replace(comparison,'\\|', ' vs ')) %>%
    filter(!is.na(p))
colors <- RColorBrewer::brewer.pal(9,'Dark2')
colors <- c(colors,'darkgrey','#D8BA90', '#7BCBD5', '#87CDA9')
pv_p <- ggplot(pval_df, aes(x=comparison, y = -log10(p), color = comparison)) + 
    geom_jitter() +
    geom_hline(yintercept = -log10(0.05), color='red') +
    labs(x = ' ', y = '-log10(Wilcox test p-value)', color = ' ')  +
    facet_wrap(~type, scale='free_x') +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    scale_color_manual(values = colors) +
    panel_border()
figurename <- file.path(figurepath, 'count_pvalue_plot.pdf')
ggsave(pv_p, file=figurename)
message('Plotted: ', figurename)


#make color
gene_count_p <- ggplot(data = plot_df,
        aes(x = map_type, y = gene_count, fill = type)) +
    geom_bar(stat='identity') +
    facet_grid(~samplename) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    scale_fill_manual(values = rev(colors)) +
    labs(x = 'Pipeline', y = 'Number of detected genes', fill= ' ')
#merge_df %>% group_by(map_type,samplename,type) %>% summarize(n())
#merge_df.spread <- merge_df %>% 
#    mutate(samplename = str_c(samplename, map_type)) %>% 
#    select(-map_type) %>% 
#    spread(samplename, abundance) %>%
#    select(grep('Sample_A_1|id|name|type', names(.))) %>% 
#    filter( Sample_A_1conventional == 0,Sample_A_1customized==0, Sample_A_1kallisto==0, Sample_A_1salmon>0)    
    
df <- merge_df %>%
    group_by(type, map_type, samplename) %>%
    summarize(read_count = sum(abundance)) %>%
    ungroup() %>%
    group_by(map_type, samplename )%>%
    do(data_frame(
        percentage_count = .$read_count/sum(.$read_count) * 100,
        type=.$type
    )) %>%
    ungroup %>%
    mutate(mix = get_mix(samplename)) %>%
    mutate(replicate = get_sample_number(samplename)) %>%
    mutate(x_name = str_c(mix, replicate)) %>%
    mutate(type = factor(type, levels = levels(plot_df$type)))%>%
    tbl_df


dist_p <-ggplot(data = df, aes(x=x_name, y = percentage_count, fill=type)) + 
    geom_bar(stat='identity') + 
    facet_grid(~map_type) +
    labs(x = ' ', y = '% Transcripts', fill= ' ') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    scale_fill_manual(values = colors)

p<-ggdraw() +
    draw_plot(gene_count_p, 0.02, 0.6, 0.95, 0.39)+
    draw_plot(ercc_lm, 0.02, 0.3, 0.95, 0.3) +
    draw_plot(ercc_r2, 0.02, 0, 0.95, 0.3) +
    draw_plot_label(label=letters[1:3], x=rep(0,3), y=c(1,0.6,0.3), size=20)

p <- plot_grid(gene_count_p, ercc_lm, 
               ncol=1, align='v', axis='l',
               labels = letters[1:2],
               label_size = 20,
               rel_heights=c(2,1))
figurepath <- file.path(project_path,'/figures')
figurename <- file.path(figurepath , 'exploring_genes.pdf')
ggsave(p, file=figurename, width = 12,height=12)
message('Plotted: ', figurename)
    

average_df <- merge_df %>%
    filter(abundance > 0.1) %>%
    mutate(samplename = str_replace(samplename,'_[123]$','')) %>%
    group_by(map_type, id, type, samplename) %>%
    summarize(abundance = mean(log2(abundance))) %>%
    ungroup()  %>%
    #select(map_type, id, type, samplename, abundance) %>%
    tbl_df
 


colors <- RColorBrewer::brewer.pal(n=8,'Dark2')
colors <- c(colors,'darkblue', 'firebrick4','darkorchid4','darkslategray3')
gene_type_cor_dfs <- function(comp1, comp2, average_df){
    comparison_label <- str_c(comp1,' vs ',comp2)
    
    total_cor <- average_df %>%
        filter(map_type %in% c(comp1, comp2)) %>%
        spread(map_type, abundance, drop=T) %>%
        set_names(c('id','type','samplename','a','b')) %>%
        group_by(samplename) %>%
        summarize(correlation = cor(a, b, use='pairwise.complete.obs')) %>%
        mutate(type='Total RNA')
        
    average_df %>%
        filter(map_type %in% c(comp1, comp2)) %>%
        spread(map_type, abundance, drop=T) %>%
        set_names(c('id','type','samplename','a','b')) %>%
        group_by(type, samplename) %>%
        summarize(correlation = cor(a, b, use='pairwise.complete.obs')) %>%
        ungroup() %>%
        rbind(total_cor) %>%
        group_by(type) %>%
        summarize(average_cor = mean(correlation),
              stdev = sd(correlation)) %>%
        mutate(comparison = comparison_label) %>%
        ungroup() %>%
        return
}
all_comparisons <- merge_df %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',map_type),1,2)) %>%
    mutate(map_type = forcats::fct_reorder(map_type, pipeline_type)) %>%
    .$map_type %>% 
    levels %>% 
    combn(2)
cor_df <- map2_df(all_comparisons[1,], 
                  all_comparisons[2,], 
                  gene_type_cor_dfs, average_df) %>%
    mutate(label = ifelse(comparison %in% c('HISAT2+featureCounts vs TGIRT-map','kallisto vs salmon'),
                          'Within type','Across type')) %>%
    mutate(type = forcats::fct_reorder(type, average_cor, min)) %>%
    mutate(type = factor(type) %>% relevel('Total RNA'))  %>%
    mutate(type_label = ifelse(type=='Total RNA', 'Total RNA', 'Individual gene type'))
    
cor_p <- ggplot(data=cor_df %>%
           mutate(comparison = forcats::fct_reorder(comparison, -average_cor, sum)), 
        aes(fill=type_label, x=type, y = average_cor)) + 
    geom_bar(stat='identity') + 
    geom_errorbar(aes(ymin=average_cor - stdev, ymax=average_cor + stdev), width=0.5) +
    facet_wrap(~comparison) + 
    scale_fill_manual(values=c('grey','red')) +
    theme(axis.text.x = element_text(angle=75, hjust=1)) +
    labs(x = 'RNA type', y = expression(rho), fill='RNA type')+
    panel_border()
figurename <- file.path(figurepath , 'genes_correlations.pdf')
ggsave(cor_p, file = figurename, width = 13, height=8)


quantile_group <- rev(c('0-25%','26-50%','51-75%','76-100%'))
spreaded_df <-merge_df %>%
    mutate(samplename = str_replace(samplename,'_[1-3]$','')) %>%
    filter(abundance > 0.1) %>%
    group_by(map_type, samplename, id) %>%
    summarize(abundance = mean(log2(abundance))) %>%
    ungroup() %>%
    group_by(id, samplename) %>%
    do(data_frame(
        map_type=.$map_type,
        abundance = .$abundance,
        mean_abundance = mean(.$abundance)
    )) %>%
    ungroup() %>%
    spread(map_type, abundance, drop=T) %>%
    drop_na() %>%
    set_names(make.names(names(.))) %>%
    tbl_df

spreaded_df %>% write_feather(file.path(project_path, '/DEgenes/tpm_table.feather'))

group_expression_df <- spreaded_df %>%
    group_by(samplename) %>%
    mutate(pc_x = ntile(mean_abundance, 4)) %>%
    ungroup() %>%
    mutate(pct_group = sapply(pc_x, function(x) quantile_group[x])) %>%
    select(-pc_x, -mean_abundance,-id) %>%
    group_by(pct_group,samplename) %>%
    nest() %>%
    mutate(cor_matrix = map(data,cor)) %>%
    mutate(cor_matrix = map(cor_matrix, function(x) {data.frame(x) %>%
                                                    tibble::rownames_to_column()})
    ) %>%
    unnest(cor_matrix) %>%
    gather(map_type, cor_value, -pct_group:-rowname) %>%
    filter(rowname != map_type) %>%
    mutate(comparison = map2(rowname, map_type,
                             function(x,y) {paste(sort(c(x,y)), collapse=' vs ')})) %>%
    mutate(comparison = unlist(comparison)) %>%
    mutate(comparison = stri_replace_all_fixed(comparison,'HISAT2.feat','HISAT2+feat')) %>%
    mutate(comparison = stri_replace_all_fixed(comparison,'TGIRT.map','TGIRT-map')) %>%
    select(-rowname, -map_type) %>%
    distinct() %>%
    ungroup() %>%
    mutate(pct_group = factor(pct_group, levels = rev(quantile_group))) %>%
    mutate(comparison = str_replace(comparison,'\\.','+')) %>%
    mutate(samplename = str_replace(samplename,'_',' ')) %>%
    group_by(pct_group, comparison) %>%
    summarize(cor_value = mean(cor_value)) %>%
    ungroup() %>%
    mutate(comp_type = ifelse(comparison %in% c('HISAT2+featureCounts vs TGIRT-map', 'kallisto vs salmon'),
                              1,2)) %>%
    mutate(comparison = forcats::fct_reorder(comparison, comp_type)) %>%
    tbl_df

expression_cor_line_plot <- ggplot(group_expression_df,
                        aes(x=pct_group, y= cor_value, 
                            color = comparison, group=comparison)) +
    geom_line(size=2, alpha=0.7) +
    labs(x = 'Expression level (top %)',y="Pearson's correlation of\nestimated expression level",
         color = ' ') +
    scale_color_manual(values = c('#D68C95', '#C79770', '#7AAD6E', '#30B1B1', '#949DD5', '#C68DC8'))
figurename <- str_c(figurepath, '/cor_expressino_line_plot.pdf')
ggsave(expression_cor_line_plot, file = figurename, width=7, height=7)
message('Plotted: ', figurename)


gene_length_cor <- spreaded_df %>%
    inner_join(gene_length_df) %>%
    mutate(gene_length_group = ntile(gene_length,4)) %>%
    select(-id,-mean_abundance,-gene_length) %>%
    group_by(gene_length_group, samplename) %>%
    nest() %>%
    mutate(cor_matrix = map(data,cor)) %>%
    mutate(cor_matrix = map(cor_matrix, function(x) {data.frame(x) %>%
                                                    tibble::rownames_to_column()})
    ) %>%
    unnest(cor_matrix) %>%
    gather(map_type, cor_value, -gene_length_group:-rowname) %>%
    filter(rowname != map_type) %>%
    mutate(comparison = map2(rowname, map_type,
                             function(x,y) {paste(sort(c(x,y)), collapse=' vs ')})) %>%
    mutate(comparison = unlist(comparison)) %>%
    select(-rowname, -map_type) %>%
    distinct() %>%
    ungroup() %>%
    mutate(gene_length_group = sapply(gene_length_group, function(x) quantile_group[5-x])) %>%
    mutate(gene_length_group = factor(gene_length_group, levels = rev(quantile_group))) %>%
    mutate(comparison = stri_replace_all_fixed(comparison,'HISAT2.feat','HISAT2+feat')) %>%
    mutate(comparison = stri_replace_all_fixed(comparison,'TGIRT.map','TGIRT-map')) %>%
    mutate(comparison = str_replace(comparison,'\\.','+')) %>%
    mutate(samplename = str_replace(samplename,'_',' ')) %>%
    group_by(comparison, gene_length_group) %>%
    summarize(cor_value = mean(cor_value)) %>%
    ungroup() %>%
    mutate(comp_type = ifelse(comparison %in% c('HISAT2+featureCounts vs TGIRT-map', 'kallisto vs salmon'),
                              1,2)) %>%
    mutate(comparison = forcats::fct_reorder(comparison, comp_type)) %>%
    tbl_df

length_cor_line_plot <- ggplot(gene_length_cor,
                        aes(x=gene_length_group, y= cor_value,
                            color = comparison, group=comparison)) +
    geom_line(size=2, alpha=0.7) +
    labs(x = 'Gene length (Long to short %)',y="Pearson's correlation of\nestimated expression level",
         color = ' ', linetype=' ') +
    scale_color_manual(values = c('#D68C95', '#C79770', '#7AAD6E', '#30B1B1', '#949DD5', '#C68DC8'))
figurename <- str_c(figurepath, '/cor_gene_length_line_plot.pdf')
ggsave(length_cor_line_plot, file = figurename, width=7, height=7)
message('Plotted: ', figurename)

cor_lines <- plot_grid(length_cor_line_plot, expression_cor_line_plot,
          ncol=1, labels = letters[1:2], label_size=20)
figurename <- str_c(figurepath, '/cor_line_plots.pdf')
ggsave(cor_lines, file = figurename, width=8, height=10)
message('Plotted: ', figurename)

venn_Df <- merge_df %>% 
    mutate(samplename = str_replace(samplename,'_[1-3]$','')) %>%
    group_by(samplename, id, name, type, map_type) %>%
    summarize(abundance = mean(abundance)) %>%
    filter(abundance > 0.1)  %>% 
    ungroup() %>%
    group_by(samplename,id,name,type) %>% 
    summarize(
              maps = str_c(unique(map_type), collapse=','),
              abundance = str_c(abundance, collapse=',')
    ) %>% 
    ungroup()%>%
    mutate(samplename = str_replace(samplename,'_',' ')) %>%
    ungroup %>% 
    tbl_df 
merge_df %>% write_feather(file.path(project_path,'/DEgenes/all_tpm_table.feather'))

colors <- RColorBrewer::brewer.pal(n=8,'Dark2')
colors <- c(colors,'darkblue', 'firebrick4','darkorchid4','darkslategray3')
sep_p <- venn_Df %>% 
    group_by(maps,type,samplename) %>% 
    summarize(count = n()) %>% 
    ungroup() %>% 
    filter(!grepl(',',maps)) %>%
    mutate(pipeline_type = ifelse(grepl('HISAT|TGIR',maps),1,2)) %>%
    mutate(maps = forcats::fct_reorder(maps, pipeline_type)) %>%
    ggplot(data =., aes(x = maps, fill=type, y = count)) + 
        geom_bar(stat='identity')  +
        facet_grid(.~samplename, scale='free_x') +
        theme(axis.text.x = element_text(angle=90, hjust= 1, vjust=0.5)) +
        labs(x =' ', fill=' ', y = 'Number of genes')+
        scale_fill_manual(values = colors)
figurename <- str_c(figurepath,'/pipeline_specific_genes.pdf')
ggsave(sep_p, file = figurename, height = 7, width=14)
message('Plotted: ', figurename)


venn_Df<-venn_Df %>%   
    group_by(maps, samplename) %>% 
    summarize(number_of_genes = n()) %>%
    ungroup()

gene_count_df <- merge_df %>% 
    filter(abundance > 0.1) %>%
    group_by(map_type, samplename) %>%
    summarize(gene_count = n()) 

fm_count <- gene_count_df %>% friedman.test(gene_count~map_type|samplename,data=.)

get_all_p <-function(compare, gene_count_df){
    comp1 <- str_split(compare,'\\|')[[1]][1]
    comp2 <- str_split(compare,'\\|')[[1]][2]
    wilcox_p <- gene_count_df  %>%
        filter(grepl(comp1,map_type, fixed=T) | grepl(comp2,map_type, fixed=T)) %>%
        arrange(samplename) %>%
        wilcox.test(gene_count~map_type, data=., paired=T) %>%
        .$p.value
    return(data.frame(wilcox_p, compare))
}

merge <- function(x,y){
    comps <- c(as.character(x),as.character(y))
    comp <- sort(comps)
    return(str_c(comp, collapse='|') )
}

pval_df <- all_comparison %>%
    data.frame() %>%
    mutate(comparison = map2(X1,X2, merge)) %>%
    select(comparison) %>%
    unnest(comparison) %>%
    distinct() %>%
    mutate(wilcox_p = map(comparison, get_all_p, gene_count_df)) %>%
    unnest(wilcox_p) %>%
    tbl_df


spreaded_merge <- merge_df %>% 
    mutate(samplename = str_replace(samplename,'_[1-3]$','')) %>%
    group_by(samplename, id, name, type, map_type) %>%
    summarize(abundance = mean(abundance)) %>%
    filter(abundance > 0.1)  %>% 
    ungroup() %>%
    mutate(abundance = 1) %>% 
    spread(map_type, abundance, fill=0) %>%
    select(-3,-4) %>%
    tbl_df

color_i <- 0
for (sm in spreaded_merge$samplename %>% unique){
    color_i <- color_i + 1
    figurename <- str_c(figurepath, '/upset_',sm,'.pdf')
    pdf(figurename, width=10,height=7,onefile=FALSE)
    upset(spreaded_merge %>% 
          filter(samplename==sm) %>% 
          select(-samplename) %>%
          data.frame,
        order.by = "freq",
        main.bar.color = colors[color_i],
        decreasing=F,
        empty.intersections = "on",
        nsets = 6, number.angles = 30, 
        point.size = 3.5, line.size = 2,  
        #c(intersection size title, intersection size tick labels, 
        #  set size title, set size tick labels, set names, numbers above bars).
        text.scale=c(1.8, 1.8, 1, 1, 1.3, 1.8),
        mb.ratio = c(0.7,0.3))
    dev.off()
    message('Plotted: ', figurename)
}



gene_count_table <- str_c(figurepath, '/gene_count.csv')
gene_count_df %>% write_csv(gene_count_table)
