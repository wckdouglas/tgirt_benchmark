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


source('ercc_correlation.R')
counts_to_tpm <- function(counts, lengths) {
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

gene_file <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/all_genes.tsv' %>%
    read_tsv() %>%
    select(gene_id, name, type) %>%
    dplyr::rename(id = gene_id) %>%
    unique()

gene_length_df <- '/stor/work/Lambowitz/ref/benchmarking/human_transcriptome/genes.length' %>%
    read_tsv()  %>%
    group_by(id) %>%
    summarize(gene_length = min(gene_length)) %>%
    ungroup() %>%
    tbl_df
 
#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
alignment_free <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
    .[!grepl('genome',.)] %>%
    .[!grepl('_[0-9]+',.)] %>%
    map_df(read_feather)  %>%
    gather(samplename, abundance, -id, -map_type) %>%
    tbl_df

# read genome mapping count files
files <- c(file.path(project_path, 'genome_mapping/Trim/conventional/counts/feature_counts.tsv'),
           file.path(project_path,'genome_mapping/Counts/RAW/combined_gene_count.tsv'))
labels <- c('conventional','customized')
genome_df <- map2(files, labels, function(x,y) read_tsv(x) %>% 
                      mutate(map_type=y) %>%
                      set_names(str_replace_all(names(.),'-','_'))) %>%
    purrr::reduce(rbind) %>%
    inner_join(gene_length_df) %>%
    mutate(id = str_replace(id,'\\-$',''))%>%
    gather(samplename, abundance, -id,-map_type, -gene_length) %>%
    mutate(id = ifelse(!grepl('ERCC',id),str_replace(id, '\\-[0-9]+$',''),id)) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]+$','') ,id)) %>%
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
                                grepl('conventional',.$map_type) ~ "HISAT2+FeatureCounts",
                                grepl('customized', .$map_type) ~ "TGIRT-map",
                                TRUE~ .$map_type)
    ) %>%
    tbl_df

genome_df %>% 
    spread(samplename,abundance) %>%
    write_feather('/stor/work/Lambowitz/cdw2854/bench_marking/DEgenes/genome_abundance_tpm.feather')

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
    mutate(samplename = str_replace(samplename, 'Sample_','')) %>%
    mutate(samplename = str_replace(samplename, '_','')) %>%
    mutate(type = forcats::fct_reorder(type, gene_count , sum))

friedman_p <-  plot_df %>% 
    group_by(map_type, samplename) %>% 
    summarize(gc = sum(gene_count)) %>% 
    friedman.test(gc~map_type|samplename, data=.) %>%
    .$p.value

friedman_per_type_p <- plot_df %>%
    group_by(type) %>%
    nest() %>%
    mutate(friedman = map(data, ~friedman.test(gene_count~map_type|samplename, data=.))) %>%
    mutate(p = map_dbl(friedman, function(x) x$p.value)) %>%
    unnest(p)

#per type p value
all_comparison <- gtools::permutations(n=4,r=2,
                                      v=unique(plot_df$map_type),
                                      repeats.allowed=F)
get_all_p <-function(x1, x2, all_comparison, plot_df){
    compare <- str_c(x1, x2, sep='|')

    wilcox_p <- plot_df %>%
        filter(map_type == x1 | map_type == x2) %>%
        group_by(map_type, samplename) %>%
        summarize(gc = sum(gene_count)) %>% 
        wilcox.test(gc~map_type, data=., paired=T) %>%
        .$p.value
    return(data.frame(wilcox_p, compare))
}
pval_df <- all_comparison %>%
    data.frame() %>%
    mutate(wilcox_p = map2(X1, X2, get_all_p, all_comparison, plot_df)) %>%
    unnest(wilcox_p) %>%
    tbl_df

#per type p value
all_comparison = gtools::permutations(n=4,r=2,
                                      v=unique(plot_df$map_type),
                                      repeats.allowed=F)
get_p <-function(i, all_comparison, plot_df){
    compare <- str_c(all_comparison[i,1],all_comparison[i,2],sep='|')

    wilcox_per_type_p <- plot_df %>%
        filter(map_type == all_comparison[i,1] | map_type == all_comparison[i,2]) %>%
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
colors <- RColorBrewer::brewer.pal(12,'Paired')
colors <- c(colors,'darkgrey')
pv_p <- ggplot(pval_df, aes(x=comparison, y = -log10(p), color = type)) + 
    geom_jitter() +
    geom_hline(yintercept = -log10(0.05), color='red') +
    labs(x = ' ', y = '-log10(Wilcox test p-value)', color = ' ')  +
    facet_wrap(~comparison, scale='free_x') +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    scale_color_manual(values = colors) +
    panel_border()
figurename <- file.path(figurepath, 'count_pvalue_plot.pdf')
ggsave(pv_p, file=figurename)
message('Plotted: ', figurename)


#make color
gene_count_p <- ggplot(data = plot_df,
        aes(x = samplename, y = gene_count, fill = type)) +
    geom_bar(stat='identity') +
    facet_grid(~map_type) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    scale_fill_manual(values = colors) +
    labs(x = ' ', y = 'Number of detected genes', fill= ' ')
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

p<-plot_grid(gene_count_p, ercc_lm, ercc_r2, 
             labels=letters[1:3], ncol=1, label_size=20)
figurepath <- '/stor/work/Lambowitz/cdw2854/bench_marking/figures'
figurename <- file.path(figurepath , 'exploring_genes.pdf')
ggsave(p, file=figurename, width = 10,height=12)
message('Plotted: ', figurename)
    



colors <- RColorBrewer::brewer.pal(n=8,'Dark2')
colors <- c(colors,'darkblue', 'firebrick4','darkorchid4','darkslategray3')
#colors <- RColorBrewer::brewer.pal(12,'Paired')
#colors <- c(colors,'darkgrey')
samples <- merge_df %>% .$samplename %>% str_replace(.,'_[123]','') %>% unique 
for (sample_type in samples){
    figurename <- str_c(figurepath, '/pair_type',sample_type,'.png')
    png(figurename,width = 18, height = 16, unit='in', res=200)
    matrix_df <- merge_df %>% 
        filter(grepl(sample_type,samplename)) %>% 
        filter(abundance > 0.1) %>% 
        group_by(map_type, id, type) %>%
        summarize(abundance = mean(log2(abundance))) %>%
        ungroup()  %>%
        spread(map_type, abundance, drop=T) %>%
        inner_join(tbl_df(.) %>% 
                   group_by(type) %>% 
                   summarize(count=n())
        ) %>%
        mutate(type = str_c(type, ' (n=', count,')')) %>%
        select(-id, -count)    %>%
        set_names(make.names(names(.)))
    p <- ggpairs(matrix_df, 
                columns = names(matrix_df%>%select(-type)),
                aes(alpha=0.7,color = type)) 
    for (i in 1:p$nrow){
        for (j in 1:p$ncol){
            if (i!=j){
                p[i,j] = p[i,j] + scale_color_manual(values=colors)
            }else{
                p[i,j] = p[i,j] + scale_fill_manual(values=colors)
            }
        }
    }
    print(p)
    dev.off()
    message('Plotted: ', figurename)
}

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

spreaded_df %>% write_feather('/stor/work/Lambowitz/cdw2854/bench_marking/DEgenes/tpm_table.feather')

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
    tbl_df

expression_cor_line_plot <- ggplot(group_expression_df, 
                        aes(x=pct_group, y= cor_value, linetype=samplename,
                            color = comparison, group=interaction(comparison,samplename))) + 
    geom_line() + 
    labs(x = 'Expression level (top %)',y="Pearson's correlation of\nestimated expression level", 
         color = ' ', linetype=' ') +
    scale_color_manual(values = RColorBrewer::brewer.pal(6,'Dark2'))
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
    tbl_df

length_cor_line_plot <- ggplot(gene_length_cor, 
                        aes(x=gene_length_group, y= cor_value, linetype=samplename,
                            color = comparison, group=interaction(comparison,samplename))) + 
    geom_line() + 
    labs(x = 'Gene length (Long to short %)',y="Pearson's correlation of\nestimated expression level", 
         color = ' ', linetype=' ') +
    scale_color_manual(values = RColorBrewer::brewer.pal(6,'Dark2'))
figurename <- str_c(figurepath, '/cor_gene_length_line_plot.pdf')
ggsave(length_cor_line_plot, file = figurename, width=7, height=7)
message('Plotted: ', figurename)

cor_lines <- plot_grid(length_cor_line_plot, expression_cor_line_plot, 
          ncol=1, labels = letters[1:2], label_size=20) 
figurename <- str_c(figurepath, '/cor_line_plots.pdf')
ggsave(cor_lines, file = figurename, width=8, height=10)
message('Plotted: ', figurename)


plot_venn <- function(which_sample,d){
    colors <- gg_color_hue(4)
    plot_df <- d %>% filter(samplename == which_sample)
    samples <- plot_df %>% filter(!grepl(',',maps)) %>% .$maps %>%unique
    sample_1 <- samples[1]
    sample_2 <- samples[2]
    sample_3 <- samples[3]
    sample_4 <- samples[4]

    figurename <- str_c(figurepath,'/venn_gene_',which_sample,'.pdf')
    figurename <- str_replace_all(figurename,' ','_')
    pdf(figurename, width = 15, height=10)
    venn.plot <- draw.quad.venn(
                area1 = filter(plot_df, grepl(sample_1,maps, fixed=T))$number_of_genes %>% sum,
                area2 = filter(plot_df, grepl(sample_2,maps, fixed=T))$number_of_genes %>% sum,
                area3 = filter(plot_df, grepl(sample_3,maps, fixed=T))$number_of_genes %>% sum,
                area4 = filter(plot_df, grepl(sample_4,maps, fixed=T))$number_of_genes %>% sum,

                n12 = filter(plot_df, grepl(sample_1,maps, fixed=T), grepl(sample_2,maps, fixed=T))$number_of_genes %>% sum,
                n13 = filter(plot_df, grepl(sample_1,maps, fixed=T), grepl(sample_3,maps, fixed=T))$number_of_genes %>% sum,
                n14 = filter(plot_df, grepl(sample_1,maps, fixed=T), grepl(sample_4,maps, fixed=T))$number_of_genes %>% sum,
                n23 = filter(plot_df, grepl(sample_2,maps, fixed=T), grepl(sample_3,maps, fixed=T))$number_of_genes %>% sum,
                n24 = filter(plot_df, grepl(sample_2,maps, fixed=T), grepl(sample_4,maps, fixed=T))$number_of_genes %>% sum,
                n34 = filter(plot_df, grepl(sample_3,maps, fixed=T), grepl(sample_4,maps, fixed=T))$number_of_genes %>% sum,

                n123 = filter(plot_df, grepl(sample_1,maps, fixed=T),
                              grepl(sample_2,maps, fixed=T), grepl(sample_3,maps, fixed=T))$number_of_genes %>% sum,
                n124 = filter(plot_df, grepl(sample_1,maps, fixed=T),
                              grepl(sample_2,maps, fixed=T), grepl(sample_4, maps, fixed=T))$number_of_genes %>% sum,
                n134 = filter(plot_df, grepl(sample_1,maps, fixed=T),
                              grepl(sample_3,maps, fixed=T), grepl(sample_4, maps, fixed=T))$number_of_genes %>% sum,
                n234 = filter(plot_df, grepl(sample_2,maps, fixed=T),
                              grepl(sample_3,maps, fixed=T), grepl(sample_4, maps, fixed=T))$number_of_genes %>% sum,

                n1234 = filter(plot_df, maps == str_c(sample_1, sample_2,
                                                      sample_3, sample_4,sep=','))$number_of_genes %>% sum,
                category = samples,
                fill = colors,
                lty = "dashed",
                cex = 1.5, cat.cex = 1.5,
                cat.col = colors)
    dev.off()
    message('Plotted: ', figurename)
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


venn_Df <- merge_df %>% 
    mutate(samplename = str_replace(samplename,'_[1-3]$','')) %>%
    group_by(samplename, id, name, type, map_type) %>%
    summarize(abundance = sum(abundance)) %>%
    filter(abundance > 0)  %>% 
    ungroup() %>%
    group_by(samplename,id,name,type) %>% 
    summarize(
              maps = str_c(unique(map_type), collapse=','),
              abundance = str_c(abundance, collapse=',')
    ) %>% 
    ungroup %>% 
    tbl_df 

sep_p <- venn_Df %>% 
    group_by(maps,type,samplename) %>% 
    summarize(count = n()) %>% 
    ungroup() %>% 
#    group_by(maps,samplename) %>% 
#    do(data_frame(
#            type = .$type, 
#            count = .$count/sum(.$count)*100
#            )
#    ) %>%  
    filter(!grepl(',',maps)) %>%
    ggplot(data =., aes(x = samplename, fill=type, y = count)) + 
        geom_bar(stat='identity')  +
        facet_grid(.~maps, scale='free_x') +
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

ps <- lapply(venn_Df$samplename%>%unique,plot_venn, venn_Df)

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


