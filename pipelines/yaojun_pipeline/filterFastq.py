#!/usr/bin/env python

import argparse
import pyximport
pyximport.install()
import filterFastqFunc

def getopt():
    help_msg = 'This is a prgram to remove/filter records from id list output is paired end fastq files'
    parser = argparse.ArgumentParser(description=help_msg)
    parser.add_argument('-1','--fq1',help='Fastq file 1',required=True)
    parser.add_argument('-2','--fq2',help='Fastq file 2', required=True)
    parser.add_argument('-i','--idfile', help='id file (one line per id)', required=True)
    parser.add_argument('-v',action='store_true')
    parser.add_argument('-o','--outprefix',required=True, help ='output prefix of the filtered fastqs')
    return parser.parse_args()

def readID(id_file):
    with open(id_file,'r') as id_f:
        id_list = {line.strip():1 for line in id_f}
    print 'Read %i unique ids' %(len(id_list))
    return id_list

def main(args):
    fastq1 = args.fq1
    fastq2 = args.fq2
    id_file = args.idfile
    inverse = args.v
    outputprefix = args.outprefix

    id_list = readID(id_file)
    seq_count = filterFastqFunc.filterFastq(fastq1, fastq2, outputprefix, id_list, inverse)
    print 'Output %i sequences' %seq_count
    return 0


if __name__ == '__main__':
    args = getopt()
    main(args)
