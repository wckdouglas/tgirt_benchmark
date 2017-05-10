#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cowplot)
library(stringr)
library(tidyr)
library(feather)
library(purrr)
library(VennDiagram)


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
 
#read alignment free abundance file from tximport
project_path <- '/stor/work/Lambowitz/cdw2854/bench_marking'
alignment_free <- project_path %>%
    file.path('DEgenes') %>%
    list.files(path = ., pattern='abundance', full.names=T) %>%
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
#    inner_join(gene_length_df) %>%
    mutate(id = str_replace(id,'\\-$',''))%>%
    gather(samplename, abundance, -id,-map_type)%>%#, -gene_length) %>%
    mutate(id = ifelse(!grepl('ERCC',id),str_replace(id, '\\-[0-9]+$',''),id)) %>%
    mutate(id = ifelse(grepl('^TR|NM|MT',id), str_replace(id,'[0-9]+$','') ,id)) %>%
    mutate(id = ifelse(id %in% c('MT-TL','MT-TS'), str_c(id,'1'), id)) %>%
    group_by(id, samplename, map_type) %>%
    summarize(
        abundance = sum(abundance)#, 
#       gene_length = mean(gene_length)
    ) %>%
    ungroup() %>%
    tbl_df

merge_df <- rbind(genome_df, alignment_free) %>%
    mutate(id = str_replace(id, '[0-9]-$','')) %>%
    inner_join(gene_file) %>%
    group_by(map_type,samplename, id, name,type) %>%
    summarize(abundance = sum(abundance)) %>%
    tbl_df


plot_df <- merge_df %>% 
    filter(abundance > 0) %>%
    group_by(type,samplename, map_type) %>%
    summarize(gene_count = n()) %>%
    ungroup() %>%
    mutate(samplename = str_replace(samplename, 'Sample_','')) %>%
    mutate(samplename = str_replace(samplename, '_','')) %>%
    mutate(type = forcats::fct_reorder(type, gene_count , sum))
#make color
colors <- RColorBrewer::brewer.pal(12,'Set3')
gene_count_p <- ggplot(data = plot_df,
        aes(x = samplename, y = gene_count, fill = type)) +
    geom_bar(stat='identity') +
    facet_grid(~map_type) +
    theme(axis.text.x = element_text(angle = 90)) +
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
    labs(x = ' ', y = '% count', fill= ' ') +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = colors)

p<-plot_grid(gene_count_p, dist_p, align='v',
             labels=letters[1:2], ncol=1)
figurepath <- '/stor/work/Lambowitz/cdw2854/bench_marking/figures'
figurename <- file.path(figurepath , 'exploring_genes.pdf')
ggsave(p, file=figurename, width = 10,height=10)
message('Plotted: ', figurename)
    

samples <- merge_df$samplename %>% str_replace(.,'_[123]','') %>% unique
for (sample_type in samples){
    figurename <- str_c(figurepath, '/pair_type',sample_type,'.png')
    png(figurename,width = 1080, height = 1080)
    matrix_df <- merge_df %>% 
        filter(grepl(sample_type,samplename)) %>%
        group_by(map_type, id, type) %>%
        summarize(abundance = log2(mean(abundance)+1)) %>%
        ungroup()  %>%
        spread(map_type, abundance) %>%
        select(-id)  
    p <- GGally::ggpairs(matrix_df, 
                columns = names(matrix_df%>%select(-type)),aes(alpha=0.7,color = type))
    print(p)
    dev.off()
    message('Plotted: ', figurename)
}


plot_venn <- function(which_sample,d){
    colors <- gg_color_hue(4)
    plot_df <- d %>% filter(samplename == which_sample)
    samples <- plot_df %>% filter(!grepl(',',maps)) %>% .$maps
    sample_1 <- samples[1]
    sample_2 <- samples[2]
    sample_3 <- samples[3]
    sample_4 <- samples[4]

    figurename <- str_c(figurepath,'/venn_gene_',which_sample,'.pdf')
    figurename <- str_replace_all(figurename,' ','_')
    pdf(figurename, width = 15, height=10)
    venn.plot <- draw.quad.venn(
                area1 = filter(plot_df, grepl(sample_1,maps))$number_of_genes %>% sum,
                area2 = filter(plot_df, grepl(sample_2,maps))$number_of_genes %>% sum,
                area3 = filter(plot_df, grepl(sample_3,maps))$number_of_genes %>% sum,
                area4 = filter(plot_df, grepl(sample_4,maps))$number_of_genes %>% sum,

                n12 = filter(plot_df, grepl(sample_1,maps), grepl(sample_2,maps))$number_of_genes %>% sum,
                n13 = filter(plot_df, grepl(sample_1,maps), grepl(sample_3,maps))$number_of_genes %>% sum,
                n14 = filter(plot_df, grepl(sample_1,maps), grepl(sample_4,maps))$number_of_genes %>% sum,
                n23 = filter(plot_df, grepl(sample_2,maps), grepl(sample_3,maps))$number_of_genes %>% sum,
                n24 = filter(plot_df, grepl(sample_2,maps), grepl(sample_4,maps))$number_of_genes %>% sum,
                n34 = filter(plot_df, grepl(sample_3,maps), grepl(sample_4,maps))$number_of_genes %>% sum,

                n123 = filter(plot_df, grepl(sample_1,maps),
                              grepl(sample_2,maps), grepl(sample_3,maps))$number_of_genes %>% sum,
                n124 = filter(plot_df, grepl(sample_1,maps),
                              grepl(sample_2,maps), grepl(sample_4, maps))$number_of_genes %>% sum,
                n134 = filter(plot_df, grepl(sample_1,maps),
                              grepl(sample_3,maps), grepl(sample_4, maps))$number_of_genes %>% sum,
                n234 = filter(plot_df, grepl(sample_2,maps),
                              grepl(sample_3,maps), grepl(sample_4, maps))$number_of_genes %>% sum,

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
    summarize(maps = str_c(unique(map_type), collapse=',')) %>% 
    ungroup %>% 
    group_by(maps, samplename) %>% 
    summarize(number_of_genes = n()) %>%
    ungroup()

ps <- lapply(venn_Df$samplename%>%unique,plot_venn, venn_Df)

