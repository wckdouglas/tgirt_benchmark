#!/usr/bin/env python

import pysam
import argparse
import parse_bam
import re
import numpy as np

def getopt():
    parser = argparse.ArgumentParser(description = 'Splitting bam file to uniquely mapped and multiple mapped (paired end only)')
    parser.add_argument('-t','--program', choices=['bowtie2','hisat2'], default='bowtie2',
            help = 'PRogram generating the bam file (default: bowtie2)')
    parser.add_argument('-i', '--inBam', default='-', help = 'bam file name, or stdin (default: stdin)' )
    parser.add_argument('-o','--outprefix', required=True, help='output prefix ($OUTPUTPREFIX.unique.bam, $OUTPUTPREFIX.multi.bam)')
    return parser.parse_args()

def softClipSize(aln):
    # compare softclip size, output maximum softclipped base on either side
    cigar = aln.cigarstring
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    clipped = cigarNum[cigarStr == 'S']
    return np.max(clipped) if len(clipped) > 0 else 0

def bowtie2split(bam, uniquebam, multibam, id_file):
    i = 0
    while True:
        try:
            read1 = bam.next()
            read2 = bam.next()
            i += 2
            assert read1.qname == read2.qname, 'Not paired end'
            if not read1.is_unmapped and not read2.is_unmapped and \
                    softClipSize(read1) < 20 and softClipSize(read2) < 20:
                if read1.mapq == 255 and read2.mapq == 255:
                    uniquebam.write(read1)
                    uniquebam.write(read2)
                else:
                    multibam.write(read1)
                    multibam.write(read3)
            else:
                id_file.write(read1.qname + '\n')
            if i % 10000000 == 0:
                print 'Parsed %i alignments' %i
        except StopIteration:
            break
    return 0

def hisat2split(bam, uniquebam, multibam, id_file):
    i = 0
    while True:
        try:
            read1 = bam.next()
            read2 = bam.next()
            i += 2
            assert read1.qname == read2.qname, 'Not paired end'
            if not read1.is_unmapped and not read2.is_unmapped and \
                    softClipSize(read1) < 20 and softClipSize(read2) < 20:
                if read1.get_tag('NH') == 1 and read2.get_tag('NH') == 1:
                    uniquebam.write(read1)
                    uniquebam.write(read2)
                else:
                    multibam.write(read1)
                    multibam.write(read2)
            else:
                id_file.write(read1.qname + '\n')
            if i % 10000000 == 0:
                print 'Parsed %i alignments' %i
        except StopIteration:
            break
    return 0

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
                bowtie2split(inbam, uniquebam, multibam, id_file)
            elif args.program == 'hisat2':
                hisat2split(inbam, uniquebam, multibam, id_file)
    return 0

if __name__ == '__main__':
    main()
