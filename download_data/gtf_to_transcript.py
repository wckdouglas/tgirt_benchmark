#!/usr/bin/env python

from gtf_to_bed import parse_extra_fields
import sys
from operator import itemgetter
import re

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
            print(line_template.format(gene_id = info_dict['gene_id'],
                                name = info_dict['gene_name'],
                                type = info_dict['gene_biotype'],
                                t_id = info_dict['transcript_id']))