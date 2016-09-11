#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import izip
import gzip

def filterFastq(fastq1, fastq2, outputprefix, id_list, inverse):
    new_fq1 = outputprefix + '_1P.fastq.gz'
    new_fq2 = outputprefix + '_2P.fastq.gz'
    cdef:
        int seq_count = 0
        str id1, seq1, qual1
        str id2, seq2, qual2
    with gzip.open(fastq1,'rb') as fq1, gzip.open(fastq2,'rb') as fq2, \
            gzip.open(new_fq1,'wb') as out_fq1, gzip.open(new_fq2, 'wb') as out_fq2:
        for (id1, seq1, qual1), (id2, seq2, qual2) in izip(FastqGeneralIterator(fq1), FastqGeneralIterator(fq2)):
            item = id_list.get(id1.split(' ')[0], None)
            if (item is not None and not inverse) or \
                    (item is None and inverse):
                out_fq1.write('@%s\n%s\n+\n%s' %(id1, seq1, qual1))
                out_fq2.write('@%s\n%s\n+\n%s' %(id2, seq2, qual2))
                seq_count += 1
    return seq_count
