#!/use/bin/env python

from __future__ import print_function
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import pandas as pd
from bs4 import BeautifulSoup as bs
import time
import sys
import numpy as np

if len(sys.argv) !=2:
    sys.exit('[usage] python %s <ref_path>' %sys.argv[0])

amino_acid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
          'SEC': 'U', 'SUP':'SUP', 'UNDET':'UND', 'IMET': 'M',
          'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def assign_tRNA_diff_gene(score_df):
    return score_df\
        .sort_values('Locus')\
        .assign(copy_id = lambda d: np.arange(d.shape[0]) +1)


def assign_tRNA_number(anticodon_df):
    new_d = anticodon_df\
        .sort_values('General tRNA Model Score') \
        .assign(tRNA_number = lambda d: d['General tRNA Model Score'].astype('category').cat.codes[::-1])\
        .assign(tRNA_number = lambda d: (d.tRNA_number - d.tRNA_number.max()) * -1 + 1)  \
        .groupby('tRNA_number',as_index=False) \
        .apply(assign_tRNA_diff_gene) \
        .assign(amino_acid = lambda d: d['Best Isotype Model'].map(lambda x: amino_acid[x.upper()])) \
        .assign(tRNA_id = lambda d: 'TR' + d.amino_acid + '-' +\
                                    d.Anticodon + \
                                    d.tRNA_number.astype('str') + '-' +\
                                    d.copy_id.astype('str')) \
        .pipe(lambda d: d[['GtRNAdb Gene Symbol','tRNAscan-SE ID','Locus','General tRNA Model Score', 'tRNA_id']]) 
    return new_d

def add_nm_tRNA(id, name):
    id = id
    if 'nm' in name:
        return 'NM' + id 

def get_table(url, name):
# using browser to get to the page
    ref_path = sys.argv[1]
    tab_name =  ref_path+'/%s_tRNA.info' %(name)


    browser = webdriver.PhantomJS()
    print('Opened browser')
    browser.get(url)
    print('Accessed webpage')
    
    # show all record
    select = Select(browser.find_element_by_name('tRNAlist_length'))   
    select.select_by_visible_text('All')
    print('Showed all')
    
    
    # get source and parse
    time.sleep(5) # wait for the webpage to load
    content = browser.page_source
    soup = bs(content, "lxml")
    html_tab = soup.find('table', id = 'tRNAlist')
    print('Scraped table')
    
    # make table
    df = pd.read_html(str(html_tab))[0] 

    df = df \
        .assign(Anticodon = lambda d: d.Anticodon.str.replace('?','N')) \
        .groupby(['Anticodon'], as_index=False)\
        .apply(assign_tRNA_number) \
        .reset_index() \
        .drop(['level_0','level_1','level_2'],axis = 1) \
        .assign(chrom = lambda d: d.Locus.str.split(':',expand=True).iloc[:,0]) \
        .assign(start = lambda d: d.Locus.str.split(':',expand=True).iloc[:,1]\
                                        .str.split('-',expand=True).iloc[:,0]) \
        .assign(end = lambda d: d.Locus.str.split('-',expand=True).iloc[:,1]\
                                        .str.split(' ',expand=True).iloc[:,0]) \
        .assign(strand = lambda d: d.Locus.str.extract('\(([+-])\)', expand=False)) \
        .assign(chrom = lambda d: d.chrom.str.replace('chr','')) \
        .rename(columns = {'General tRNA Model Score':'score'}) \
        .assign(tRNA_name = lambda d: d.tRNA_id.str.replace('[0-9]+-[0-9]+$','')) \
        .assign(type = 'tRNA') \
        .pipe(lambda d: d[['chrom','start', 'end', 'tRNA_name',
                           'score','strand','type','tRNA_id','GtRNAdb Gene Symbol']])

#        .assign(tRNA_id = lambda d: map(add_nm_tRNA, d.tRNA_id, d.tRNA_name))
    df.to_csv(tab_name,sep='\t',index=False)
    print('Written %s' %tab_name)
    browser.close()
    browser.quit()


url = 'http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/Hsapi38-gene-list.html'
get_table(url, 'hg38')
