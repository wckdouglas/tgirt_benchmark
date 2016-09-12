#!/usr/bin/env python

import pysam
import numpy as np
import sys
import argparse

def getopt():
    parser = argparse.ArgumentParser(description='Process multiply mapped bam and get either shortest insert size or ribosomal reads')
    parser.add_argument('-i','--infile', required=True, help = 'Input bam/sam file (use < - > for stdin)')
    parser.add_argument('-o','--outfile', required=True, help = 'Output bam/sam file (use < - > for stdout)' )
    parser.add_argument('-b','--bam_in', action='store_true', help = 'Provide this flag if bam instead of sam is used as input' )
    parser.add_argument('-z','--bam_out', action='store_true', help = 'Provide this flag if bam is needed for output')
    return parser.parse_args()

def riboChrom(name):
    return True if 'gi' in name else False

def riboChromVector(chroms):
    return np.array(map(riboChrom, chroms))

def filterAlingments(read1, read2):
    chroms = []
    isizes = []
    read1_group = []
    read2_group = []
    for r1, r2 in zip(read1, read2):
        if r1.reference_name == r2.reference_name and abs(r1.isize) == abs(r2.isize):
            read1_group.append(r1)
            read2_group.append(r2)
            isizes.append(abs(r1.isize))
            chroms.append(r1.reference_name)
    chroms = np.array(chroms)
    isizes = np.array(isizes)

    #start filtering
    size_bool = (isizes == np.min(isizes))
    ribo_bool = riboChromVector(chroms)

    read1 , read2 = map(np.array, [read1, read2])

    if len(isizes[size_bool]) == 1:
        return read1[size_bool][0], read2[size_bool][0]
    elif len(chroms[ribo_bool]) == 1:
        return read1[ribo_bool][0], read2[ribo_bool][0]
    else:
        return read1[size_bool][0], read2[size_bool][0]

def putInGroup(read1, read2, alignment):
    if alignment.is_read2:
        read2.append(alignment)
    else:
        read1.append(alignment)
    return read1, read2

def printProgress(read_count):
    if read_count % 1000 == 0:
        sys.stderr.write('Processed: %i\n' %read_count)
    return 0

def processBam(in_bam, out_bam, bam_in_bool, bam_out_bool):
    read1 = []
    read2 = []
    read_flag = 'rb' if bam_in_bool else 'r'
    write_flag = 'wb' if bam_out_bool else 'w'
    sys.stderr.write('Start processing bam file: %s\n' %(in_bam))
    sys.stderr.write('Writing to: %s\n' %(out_bam))
    with pysam.Samfile(in_bam, read_flag) as in_sam:
        with pysam.Samfile(out_bam, write_flag, template = in_sam) as out_sam:
            read_count = 0
            while True:
                try:
                    #initiate new read groups
                    alignment = in_sam.next()
                    read_id = alignment.qname
                    read1, read2 = putInGroup(read1, read2, alignment)
                    read_count += 1
                    printProgress(read_count)

                    # search same id reads
                    while True:
                        new_align = in_sam.next()
                        read_count += 1
                        printProgress(read_count)
                        if new_align.qname != read_id:
                            read1_aln, read2_aln = filterAlingments(read1, read2)
                            out_sam.write(read1_aln)
                            out_sam.write(read2_aln)
                            read1 = []
                            read2 = []
                            read1, read2 = putInGroup(read1, read2, new_align)
                            read_id = alignment.qname
                            break
                        else:
                            read1, read2 = putInGroup(read1, read2, new_align)
                except StopIteration:
                    read1_aln, read2_aln = filterAlingments(read1, read2)
                    out_sam.write(read1_aln)
                    out_sam.write(read2_aln)
                    break
    return 0

def main():
    args = getopt()
    in_bam = args.infile
    out_bam = args.outfile
    bam_in_bool = args.bam_in
    bam_out_bool = args.bam_out
    processBam(in_bam, out_bam, bam_in_bool, bam_out_bool)
    return 0

if __name__ == '__main__':
    main()
