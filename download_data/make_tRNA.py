#!/usr/bin/env python

from Bio import Entrez
import time
import os
import re
from itertools import izip
import pandas as pd
import numpy as np
import urllib2
import cjson
import sys
from collections import defaultdict

def undetermined_tRNA():
    rna_central_id = 'URS000071AE84,URS0000757F07,URS00007581F8,URS0000756537,URS0000734F3E' + \
                ',URS000073F00E,URS0000732741,URS0000756CF3,URS000072FC7F,URS000074598E'
    gene_name = 'tRNA-Und-NNN-2-1,tRNA-Und-NNN-9-1,tRNA-Und-NNN-10-1,tRNA-Und-NNN-4-1,'+\
            'tRNA-Und-NNN-5-1,tRNA-Und-NNN-7-1,tRNA-Und-NNN-3-1,tRNA-Und-NNN-8-1'+\
            ',tRNA-Und-NNN-6-1,tRNA-Und-NNN-1-1'
    return rna_central_id.split(','), gene_name.split(',')

def download_table(tablename):
    link = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_rna_central_ids&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit'
    download_command = "curl '%s' >  %s" %(link, tablename)
    os.system(download_command)

def make_table(tablename):
    df = pd.read_table(tablename)\
            .pipe(lambda d: d[d['Approved Name']\
                              .str\
                              .contains('transfer RNA|mitochondrially encoded tRNA|transfer RNA-')])\
            .pipe(lambda d: d[['Approved Name','Approved Symbol','HGNC ID',
                               'RNAcentral ID(supplied by RNAcentral)']])\
            .pipe(lambda d: d[d['RNAcentral ID(supplied by RNAcentral)'].notnull()])
    return df


def fetch_gene(rna_central_id, name):
    fasta_record = ''
    base_url = 'http://rnacentral.org/api/v1/rna/'
    try:
        html = urllib2.urlopen(base_url + '/' + rna_central_id)
        html_file = html.read()
        rna_record = cjson.decode(html_file)
        seq = rna_record['sequence']
    except (TypeError,urllib2.HTTPError):
        print 'Wrong ID:' + str(rna_central_id)+'\t'+name
        pass
    return seq

def make_name(seq_name):
    return re.sub('-[0-9]+$','',seq_name)

def make_seq(seq_dict):
    record = []
    for seq, names in seq_dict.iteritems():
        description = ','.join(names)
        names = sorted(set(map(make_name, names)))
        seq_id = names[0]
       
        record.append((seq_id, description, seq, names, len(names)))
    return record

def download_seq(tRNA_table, fastaname):
    out_count = 0
    seq_dict = defaultdict(list)
    for i, record in tRNA_table.iterrows():
        rna_central_id = str(record['RNAcentral ID(supplied by RNAcentral)'])
        name = record['Approved Symbol']
        seq = fetch_gene(rna_central_id, name)
        seq_dict[seq].append(name)
        out_count += 1
    undetermined_id, undetermined_name = undetermined_tRNA()
    for rna_central_id, name in izip(undetermined_id, undetermined_name):
        fasta_record = fetch_gene(rna_central_id, name)
        seq_dict[seq].append(name)
        out_count += 1
    print 'fetched %i tRNA records' %(out_count)
    
    out_seq_name = set()
    records = make_seq(seq_dict)
    records = sorted(records, key= lambda x: x[4])
    with open(fastaname,'w') as out_fasta:
        for i, (seq_id, description, seq, names, count) in enumerate(records):
            seq = seq.replace('U','T')
            if seq_id not in out_seq_name:
                out_seq_name.add(seq_id)
                fasta_record = '>%s\t%s\n%s' %(seq_id, description, seq)
                out_fasta.write(fasta_record + '\n')
            else:
                try:
                    seq_id = names[1]
                    out_seq_name.add(seq_id)
                    fasta_record = '>%s\t%s\n%s' %(seq_id, description, seq)
                    out_fasta.write(fasta_record + '\n')
                except IndexError:
                    print 'Index collision %s' %seq_id
    print 'Written %i tRNA records' %(i)

def main():
    tablename = os.environ['REF'] + '/human_transcriptome/tRNA_nomenclature.tsv'
    fastaname = os.environ['REF'] + '/human_transcriptome/tRNA.fa'
    #download_table(tablename)
    tRNA_table = make_table(tablename)
    download_seq(tRNA_table,fastaname)

if __name__ == '__main__':
    main()
