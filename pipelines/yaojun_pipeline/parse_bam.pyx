import pysam
from pysam.calignedsegment cimport AlignedSegment
import re
import numpy as np

def softClipSize(aln):
    # compare softclip size, output maximum softclipped base on either side
    cdef str cigar
    cigar = aln.cigarstring
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    clipped = cigarNum[cigarStr == 'S']
    return np.max(clipped) if len(clipped) > 0 else 0

def bowtie2split(bam, uniquebam, multibam, id_file):
    cdef:
        int i = 0
        AlignedSegment read1
        AlignedSegment read2
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
                    multibam.write(read2)
            else:
                id_file.write(read1.qname + '\n')
            if i % 100000 == 0:
                print 'Parsed %i alignments' %i
        except StopIteration:
            break
    return 0

def hisat2split(bam, uniquebam, multibam, id_file):
    cdef:
        int i = 0
        AlignedSegment read1
        AlignedSegment read2
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
            if i % 100000 == 0:
                print 'Parsed %i alignments' %i
        except StopIteration:
            break
    return 0
