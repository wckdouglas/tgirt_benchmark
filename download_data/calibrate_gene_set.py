#!/usr/bin/env python


import pandas as pd
import os, sys
import re

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <ref path>' %sys.argv[0])
    
gene_path = sys.argv[1]

def transcript_table():
    tablename = gene_path + '/transcripts.tsv'
    return pd.read_table(tablename) \
            .pipe(lambda d: d[['gene_id','name','type','t_id']]) \
            .drop_duplicates()

def genes_bed():
    tablename = gene_path + '/genes.bed'
    return pd.read_table(tablename,
                         names = ['chrom','start','end','gene_name',
                                  'score','strand','gene_type','gene_id'])\
            .pipe(lambda d: d[['gene_name','gene_type','gene_id']])



ncRNA = ["sense_intronic","3prime_overlapping_ncRNA",'processed_transcript',
        'sense_overlapping','Other_lncRNA', 'macro_lncRNA','non_coding','known_ncrna',
        'lincRNA','bidirectional_promoter_lncRNA', 'ribozyme','3prime_overlapping_ncrna']
smncRNA = ['misc_RNA','snRNA','piRNA','scaRNA','sRNA','scRNA']
large_rRNA = ['28S_rRNA','18S_rRNA']
small_rRNA = ['rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA']
protein_coding = ['protein_coding','TR','IG']
def changeType(x):
    type = ''
    if x in ncRNA or re.search('pseudo',x):
        type = 'Other ncRNA'
    elif re.search('TR|IG|protein',x):
        type = 'Protein coding'
    elif x.startswith('Mt_'):
        type = 'Mt'
    elif x == 'tRNA':
        type = 'tRNA'
    elif x in small_rRNA or x in large_rRNA:
        type = 'rRNA'
    elif x in smncRNA:
        type = 'Other sncRNA'
    elif x =='antisense':
        type = 'Antisense'
    else:
        type = x
    return type
    
def union_type(x, y):
    out_type=''
    types = filter(lambda s: s != 'nan', list(set(map(str, [x,y]))))
    if len(types):
        out_type = changeType(types[0])
    else:
        out_type = ','.join(types)
    return out_type

def mt_tRNA(name,in_type):
    type = ''
    if 'MT' in str(name) and in_type =='tRNA':
        type = 'Mt'
    elif (str(name).startswith('VT') or name == 'Vault') and in_type == 'Other sncRNA':
        type = 'vaultRNA'
    else:
        type = in_type
    return type

def rRNA_name(name,type):
    name = str(name)
    if type == 'rRNA':
        if '5S' in name:
            name = '5S_rRNA'
        elif '5_8S_r' in name or '5-8S' in name:
            name = '5.8S_rRNA'
    else:
        name = name
    return name
  
def fill_name(name, gene_name):
    if name == 'nan':
        name = gene_name
    return name

t_tab = transcript_table()
bed = genes_bed()
df = t_tab.merge(bed, how='outer', on ='gene_id') \
    .assign(type = lambda d:map(union_type, d.type, d.gene_type)) \
    .assign(t_id = lambda d: d.t_id.fillna('NA')) \
    .assign(type = lambda d: map(mt_tRNA, d.name, d.type))\
    .assign(name = lambda d: map(rRNA_name, d.name, d.type)) \
    .assign(name = lambda d: map(fill_name, d.name, d.gene_name))\
    .drop(['gene_name','gene_type'], axis=1) 
union_gene = gene_path + '/all_genes.tsv'
df.to_csv(union_gene, index=False, sep='\t')
print 'Written %s' %union_gene
