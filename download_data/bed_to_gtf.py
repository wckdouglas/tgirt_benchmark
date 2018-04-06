#!/usr/bin/env python

from __future__ import print_function
import sys
from operator import itemgetter

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <bed_file>' %(sys.argv[0]))


def make_line(chrom ,start, end, strand, gene_name, gene_id, gene_type):
    out_line_gene = '{chrom}\tcustomized\tgene\t{start}\t{end}\t.\t{strand}\t.\t' \
            'gene_id "{gene_id}_gene"; gene_version "1"; gene_name "{gene_name}"; ' \
            'gene_source "customized"; gene_biotype "{gene_type}";'\
                .format(chrom=chrom, start= start, end= end, 
                        strand = strand, gene_id = gene_id, 
                        gene_name = gene_name, gene_type = gene_type)

    out_line_transcript = '{chrom}\tcustomized\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t' \
            'gene_id "{gene_id}_gene"; gene_version "1"; transcript_id "{transcript_id}"; ' \
            'transcript_version "1"; gene_name "{gene_name}"; ' \
            'gene_source "customized"; gene_biotype "{gene_type}"; transcript_name "{transcript_name}"; '\
            'transcript_biotype "{transcript_biotype}";' \
                .format(chrom=chrom, start= start, 
                        end= end, strand = strand,
                        transcript_id = gene_id ,gene_id = gene_id,
                        gene_name = gene_name,
                        gene_type = gene_type, transcript_name = gene_id,
                        transcript_biotype = gene_type) 

    if 'ERCC-' in gene_id or gene_id.startswith('TR'):
        end = int(end) - 1
    out_line_exon = '{chrom}\tcustomized\texon\t{start}\t{end}\t.\t{strand}\t.\t' \
            'gene_id "{gene_id}_gene"; gene_version "1"; transcript_id "{transcript_id}"; ' \
            'transcript_version "1"; exon_number "1"; gene_name "{gene_name}"; ' \
            'gene_source "customized"; gene_biotype "{gene_type}"; transcript_name "{transcript_name}"; '\
            'transcript_biotype "{transcript_biotype}"; ' \
            'exon_id "exon_id"; exon_version "1"; transcript_support_level "1";'\
                .format(chrom=chrom, start= start, end= end, 
                        strand = strand, transcript_id = gene_id ,
                        gene_id = gene_id, gene_name = gene_name,
                        gene_type = gene_type, transcript_name = gene_id,
                        transcript_biotype = gene_type, 
                        exon_id = gene_id + '_exon')

    return out_line_gene + '\n' + out_line_transcript + '\n' + out_line_exon


with open(sys.argv[1],'r') as in_bed:
    for line in in_bed:
        fields = line.strip().split('\t')
        chrom ,start, end, strand, gene_name, gene_id, gene_type = itemgetter(0,1,2,5,3,7,6)(fields)
        line = make_line(chrom ,start, end, strand, gene_name, gene_id, gene_type)
        print(line)

