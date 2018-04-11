#!/usr/bin/env python

from gtf_to_bed import parse_extra_fields
import sys
from operator import itemgetter
import re


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


def mt_tRNA(name,in_type):
    type = ''
    if 'MT-' in str(name) or ('MT' in name and in_type =='tRNA'):
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


#gtf = '/stor/work/Lambowitz/ref/benchmarking/GRCH38_genome/genes.gtf'
if len(sys.argv) != 2:
    sys.exit('[usage] python %s <gtf file>' %sys.argv[0])
gtf = sys.argv[1]

line_template = '{gene_id}\t{name}\t{t_id}\t{type}'
header = re.sub('\{|\}','',line_template)
print(header)
for line in open(gtf, 'r'): 
    if not line.startswith('#'):
        fields = line.strip().split('\t')
        if fields[2] == "transcript":
            extra_fields = fields[8]
            info_dict = parse_extra_fields(extra_fields)
            name = info_dict['gene_name']
            gene_type = changeType(info_dict['gene_biotype'])

            name = rRNA_name(name, gene_type)
            gene_type = mt_tRNA(name, gene_type)

            gene_type = 'rRNA' if gene_type == "rDNA" else gene_type
            gene_id = info_dict['gene_id']
            if 'rRNA' in name or gene_type == 'tRNA':
                gene_id = name

            print(line_template.format(gene_id = gene_id,
                                name =  name ,
                                type =  gene_type,
                                t_id = info_dict['transcript_id']))
            
            if gene_type == "Mt" and 'MT-T' in name:
                print(line_template.format(gene_id = name,
                    name =  name ,
                    type =  gene_type,
                    t_id = info_dict['transcript_id'])) 
