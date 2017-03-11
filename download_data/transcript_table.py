#!/usr/bin/env python

import pandas as pd
import sys

def extractFields(info_fields):
    info_dict = {}
    for field in info_fields:
        field = field.strip(' ')
        attribute = field.split(' ')[0]
        try:
            value = field.split(' ')[1].strip('"')
        except IndexError:
            print info_fields
        info_dict[attribute] = value
    return info_dict

def extractFeatures(gene_line):
    fields = gene_line.strip('\n').split('\t')
    if gene_line[0]!='#' and fields[2] == 'transcript':
        info = fields[-1]
        info_fields = info.split(';')[:-1]
        info_dict = extractFields(info_fields)
        try:
            return_tuple = (info_dict['gene_id'],info_dict['gene_name'],
                info_dict['transcript_id'],info_dict['gene_biotype'])
        except KeyError:
            return_tuple = (info_dict['gene_id'],info_dict['gene_id'],
                info_dict['transcript_id'],info_dict['gene_biotype'])
#            print info_dict
#            sys.exit('Wrongly parsed')
        return return_tuple


gene_path = sys.argv[1]
gene_file = gene_path + '/genes.gtf'
t_file = gene_path + '/transcripts.tsv'
gene_id = []
t_id = []
biotype = []
with open(gene_file,'r') as genes:
    d = map(extractFeatures, genes)
d = filter(None, d)
pd.DataFrame(d, columns = ['gene_id','name','t_id','type'])\
        .drop_duplicates() \
        .to_csv(t_file, sep='\t',index=False)
print 'Written: %s' %t_file


