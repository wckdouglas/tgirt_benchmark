#!/usr/bin/env python

import pysam
import argparse
import re
import pyximport
pyximport.install(setup_args={'include_dirs': pysam.get_include()})
import parse_bam
import numpy as np

def getopt():
    parser = argparse.ArgumentParser(description = 'Splitting bam file to uniquely mapped and multiple mapped (paired end only)')
    parser.add_argument('-t','--program', choices=['bowtie2','hisat2'], default='bowtie2',
            help = 'PRogram generating the bam file (default: bowtie2)')
    parser.add_argument('-i', '--inBam', default='-', help = 'bam file name, or stdin (default: stdin)' )
    parser.add_argument('-o','--outprefix', required=True, help='output prefix ($OUTPUTPREFIX.unique.bam, $OUTPUTPREFIX.multi.bam)')
    return parser.parse_args()

def main():
    args = getopt()
    outprefix = args.outprefix
    uniqueBam = outprefix + '.unique.bam'
    multiBam = outprefix + '.multi.bam'
    unmapped_id = outprefix + '.id.dat'
    with pysam.Samfile(args.inBam, 'rb') as inbam:
        with pysam.Samfile(uniqueBam,'wb', template=inbam) as uniquebam,  \
                pysam.Samfile(multiBam,'wb',template=inbam) as multibam, \
                open(unmapped_id,'w') as id_file:
            if args.program == 'bowtie2':
                parse_bam.bowtie2split(inbam, uniquebam, multibam, id_file)
            elif args.program == 'hisat2':
                parse_bam.hisat2split(inbam, uniquebam, multibam, id_file)
    return 0

if __name__ == '__main__':
    main()
