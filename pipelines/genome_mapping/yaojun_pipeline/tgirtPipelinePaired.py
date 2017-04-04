#!/bin/env  python
# This is the pair-end pipeline for tgirt sequencing
# Mapping with hisat + bowtie local
# and extract tRNA reads for reassigning counts


from __future__ import division
from multiprocessing import Pool
import os
import sys
import time
import glob
import argparse

def makeFolder(folder):
    """
        Input a folder name and make a folder if it is non-existed
    """
    sys.stderr.write('Creating %s....\n' %folder)
    if os.path.isdir(folder):
        sys.stderr.write('%s exists.\n' %folder)
    else:
        os.mkdir(folder)
        sys.stderr.write('Created %s.\n' %folder)
    return 0

def runProcess(args):
    command, samplename = args
    print '[%s] Running: %s' %(samplename, command)
    start = time.time()
    os.system(command)
    end = time.time() - start
    sys.stderr.write('[%s] Used time %.3f min\n' %(samplename, end/60))
    return 0

def trimmomatic(fastqfile,resultpath,sampleName,cores,adaptors):
    """ trimmomatic commands
        input is: fastqFile in gz or unzip form
        result directory
        sample name
        and cores
    """
    options='ILLUMINACLIP:%s:2:10:10:1:true ' %adaptors+ \
            'LEADING:10 TRAILING:10  SLIDINGWINDOW:4:8 ' +\
            'MINLEN:10  AVGQUAL:20'
    resultfile=resultpath+'/'+sampleName + '.fastq.gz'
    command='trimmomatic PE ' +\
            '-threads %i -phred33 -basein %s ' %(cores,fastqfile) +\
            '-baseout %s %s' %(resultfile,options)
    runProcess((command,sampleName))
    return 0

def hisat_pariedEnd(datapath,result_dir,sampleName,cores,index, spliceFile):
    """
        this run hisat in single end with default settings,
        input:
                fastq file
                result directory
                sample name
                core to use
                bowtie index
    """
    file1 = datapath + '/' + sampleName + '_1P.fastq.gz'
    file2 = datapath + '/' + sampleName + '_2P.fastq.gz'
    start = time.time()
    outFile = '%s/%s' %(result_dir, sampleName)
    command = 'hisat2 ' +\
            '--known-splicesite-infile %s --no-discordant --no-mixed -k 10 ' %(spliceFile) + \
            '--threads %i -x %s -1 %s -2 %s ' %(cores, index, file1, file2) + \
            '| samtools view -b -@ %i ' %(cores) + \
            '> %s.bam'  %outFile
    runProcess((command,sampleName))
    split_command = 'python splitBam.py --outprefix=%s --program=hisat2 --inBam=%s.bam' %(outFile,outFile)
    runProcess((split_command,sampleName))
    return outFile

def fastqRemoveID(sampleName,fastqPath,resultpath,idFile, cores):
    """
        filter fastq file using id file for pair end reads
        input:
            sample name
            input fastq path
            result fastq path
            id file
    """
    commands = ["filterFastq -v -q  %s/%s_%dP.fastq.gz -i %s > %s/%s_%dP.fastq"\
            %(fastqToolsPath,fastqPath,sampleName,end,\
            idFile,resultpath,sampleName,end) for end in [1,2]]
    Pool(cores).map(runProcess,[(command, sampleName) for command in commands])
    return 0

def bowtie_pairedEnd(datapath,resultpath,index,sampleName,cores):
    """
        run bowtie on unmapped reads after bowtie local mapping
        input:
            fastq file
            result directory
            bowtie index
            sample name
            core to use
    """
    file1 = datapath + '/' + sampleName + '_1P.fastq.gz'
    file2 = datapath + '/' + sampleName + '_2P.fastq.gz'
    resultfile = '%s/%s' %(resultpath,sampleName)
    command = 'bowtie2 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 --threads %i ' %(cores)+ \
            '--time --no-mixed --no-discordant -x %s -1 %s -2 %s ' %(index,file1,file2) + \
            '| samtools view -@ %i -b - ' %(cores) + \
            '> %s.bam' %(resultfile)
    runProcess((command,sampleName))
    split_command = 'python splitBam.py --outprefix=%s --inBam=%s.bam --program=bowtie2' %(resultfile,resultfile)
    runProcess((split_command,sampleName))
    return resultfile

def mergeBams(bam1,bam2,resultPath,sampleName,cores):
    """
        This merge two bam files and sort them
        input
            first bam file
            second bam file
            result directory
            sample name
            core to use
    """
    outBam = resultPath + '/' + sampleName + '.bam'
    command = 'samtools cat %s %s ' %(bam1,bam2)+\
	    '| samtools fixmate -O bam -r - - ' +\
        '> ' + outBam
    runProcess((command,sampleName))
    return outBam

def multibamToPrimary(hisat_multi_bam, bowtie2_multi_bam, hisat_uniq_bam,
                    bowtie2_uniq_bam, resultpath, sampleName, cores):
    "selecting best from multiply mapped bam reads and combine with uniq map reads"
    combined_bam = '%s/%s.primary.bam' %(resultpath, sampleName)
    command = 'samtools cat %s %s ' %(hisat_multi_bam, bowtie2_multi_bam) +\
        '| python reduce_multi_reads.py --infile=- --outfile=- ' +\
            '--bam_in --bam_out ' +\
        '| samtools cat %s %s - ' %(hisat_uniq_bam, bowtie2_uniq_bam) +\
        '| samtools sort -@ %i -n -O bam -T %s/%s' %(cores, resultpath, sampleName) +\
        '>  %s ' %(combined_bam)
    runProcess((command,sampleName))
    return combined_bam

def getBedMappedID(uniqueBam, idpath, sampleName, bed_ref):
    idFile = '%s/%s.id.dat' %(idpath,sampleName)
    command = 'bedtools intersect -s -abam %s -b %s -f 0.5' %(uniqueBam, bed_ref) + \
            "| samtools view " +\
            "| awk {'print $1'} > %s" %idFile
    runProcess((command,sampleName))
    return idFile


def id2Fastq(fastqPath,sampleName,idFile,resultpath, cores):
    """
        Using ID file to extract sequence from paired end fastq files
        input:
            fastq directory
            sample name
            id file
            result fastq path
    """
    commands = ["seqtk subseq %s/%s_%dP.fastq.gz %s |gzip > %s/%s_%dP.fastq.gz" \
            %(fastqPath,sampleName,end,\
            idFile,resultpath,sampleName,end) for end in [1,2]]
    Pool(cores).map(runProcess,[(command, sampleName) for command in commands])
    return 0

def rnaRemap(datapath,cores,sampleName, countpath,resultpath,rna_index,strand, rna_type):
    """
        Mapping all unmapped + tRNA reads to tRNA reference and
        extract the reads that mapped to tRNA locus
        input:
            fastq file (unmapped + tRNA)
            core to use
            sample name
            result directory
    """
    start = time.time()
    file1 = datapath + '/' + sampleName + '_1P.fastq.gz'
    file2 = datapath + '/' + sampleName + '_2P.fastq.gz'
    out_bam = resultpath + '/' + sampleName + '.%s.bam' %(rna_type)
    out_bed = resultpath + '/' + sampleName + '.%s.bed' %(rna_type)
    if strand == 0:
        type = ' --norc '
    elif strand == 1:
        type = ' --nofw '
    elif strand == 2:
        type = ' '
    map_command = 'bowtie2 --threads %i ' %(cores)+ \
            '--local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 --no-mixed --no-discordant ' +\
            type + '-x %s -1 %s -2 %s ' %(rna_index, file1,file2) + \
            "| samtools view -bF4 -@ %i "%(cores) + \
            '| tee %s ' %out_bam +\
            '| bedtools bamtobed -mate1 -bedpe -i - '+\
            '| python bedpetobed.py /dev/stdin ' +\
            '> %s' %(out_bed)
    count_command = 'cat %s' %(out_bed) +\
        '| cut -f1' + \
        '| sort ' +\
        '| uniq -c ' +\
        "| awk '{print $2,$1}' OFS='\\t'" +\
        '> %s/%s.%s.counts' %(countpath, sampleName, rna_type)
    runProcess((map_command, sampleName))
    runProcess((count_command, sampleName))
    return out_bam


def countBam(bamFile, bedpath, countpath, sampleName):
    snc_command = 'bedtools pairtobed -s -f 0.5 -abam %s '%(bamFile) +\
            '-b %s/sncRNA_no_tRNA.bed' %(bedpath) +\
        '| bedtools bamtobed -mate1 -bedpe -i - ' +\
        '| python bedpetobed.py /dev/stdin ' +\
        '| bedtools  coverage -s -counts -F 0.1 -a %s/sncRNA_no_tRNA.bed -b -' %(bedpath) +\
        '| cut -f4,7-9 ' +\
        '> %s/%s.sncRNA.counts' %(countpath, sampleName)
    non_snc_command = 'bedtools pairtobed -s -f 0.5 -type neither '+\
            '-abam %s -b %s/sncRNA_rRNA.bed ' %(bamFile, bedpath) +\
        '| bedtools bamtobed -mate1 -bedpe -i - ' +\
        '| python bedpetobed.py /dev/stdin ' +\
        '| bedtools coverage -s -counts -F 0.1 -a %s/genes_no_sncRNA_rRNA_tRNA.bed -b - ' %(bedpath) +\
        '| cut -f4,7-9 ' +\
        '> %s/%s.non_sRNA.counts' %(countpath, sampleName)
    runProcess((snc_command, sampleName))
    runProcess((non_snc_command,sampleName))
    return 0

def programSequenceControl(fastqFile,resultpath,cores,humanIndex, bedpath,
                           tRNA_index, rRNA_index, adaptors,spliceFile,strand):
    sampleName = '_'.join(fastqFile.split('/')[-1].split('_')[:-2])
    # set sample name and result path
    start = time.time()
    trimResultPath = resultpath + '/trimmed'
    hisatResultPath = resultpath + '/hisat'
    hisatUnmappedPath = hisatResultPath + '/unmapped'
    hisatMappedPath = hisatResultPath + '/mappedBam'
    bowtieResultPath = resultpath + '/bowtie2'
    bowtieMappedPath = bowtieResultPath + '/mapped'
    mergedPath = resultpath + '/mergeBam'
    mergedBamPath = mergedPath + '/bamFiles'
    uniqueBamPath = mergedBamPath + '/uniqueBam'
    multiBamPath = mergedBamPath + '/multiBam'
    mergedBedPath = mergedPath + '/bedFiles'
    uniqueBedPath = mergedBedPath + '/uniqueBed'
    multiBedPath = mergedBedPath + '/multiBed'
    countPath = mergedPath + '/countFiles'
    alltRNAfastqPath = mergedPath + '/alltRNA'
    allrRNAfastqPath = mergedPath + '/allrRNA'

    # make result directories
    folders = [resultpath,trimResultPath,hisatResultPath,hisatUnmappedPath,hisatMappedPath,
            bowtieResultPath,bowtieMappedPath,mergedPath,mergedBamPath,uniqueBamPath, multiBamPath,
            countPath, alltRNAfastqPath, mergedBedPath, uniqueBedPath, multiBedPath,
            allrRNAfastqPath]
    map(makeFolder,folders)

    # trimmomatic
    trimmomatic(fastqFile, trimResultPath,sampleName,cores,adaptors)

    # hisat map
    hisatOutprefix = hisat_pariedEnd( trimResultPath,
            hisatMappedPath,
            sampleName, cores ,
            humanIndex, spliceFile)

    # extract unmapped sequence
    id2Fastq(trimResultPath, sampleName, hisatOutprefix + '.id.dat', hisatUnmappedPath, cores)

    # local mapped
    bowtieOutprefix = bowtie_pairedEnd(hisatUnmappedPath, bowtieResultPath,
                                        humanIndex, sampleName,cores)

    #merge mapped bams
    primary_bam = multibamToPrimary(hisatOutprefix + '.multi.bam', bowtieOutprefix + '.multi.bam',
                                    hisatOutprefix + '.unique.bam', bowtieOutprefix + '.unique.bam',
                                        uniqueBamPath, sampleName,cores)

    #all tRNA ID
    tRNA_id_file = getBedMappedID(primary_bam, alltRNAfastqPath, sampleName, bedpath + '/tRNA.bed')
    rRNA_id_file = getBedMappedID(primary_bam, allrRNAfastqPath, sampleName, bedpath + '/rRNA.bed')
    id2Fastq(trimResultPath, sampleName, tRNA_id_file, alltRNAfastqPath, cores)
    id2Fastq(trimResultPath, sampleName, rRNA_id_file, allrRNAfastqPath, cores)

    #mapping tRNA
    tRNA_bam = rnaRemap(alltRNAfastqPath,cores,sampleName, countPath, alltRNAfastqPath, tRNA_index, strand, 'tRNA')
    rRNA_bam = rnaRemap(allrRNAfastqPath,cores,sampleName, countPath, allrRNAfastqPath, rRNA_index, 2 ,'rRNA')

    #countBAm
    countBam(primary_bam, bedpath, countPath, sampleName)

    #count tRNA

    usedTime = time.time()-start
    print 'Finished: %s in %.3f hr ' %(sampleName ,usedTime/3600)
    return 0

def getopt():
    parser = argparse.ArgumentParser(description='Pipeline for mapping and counting for TGIRT-seq paired end data')
    parser.add_argument('-q', '--fastq', help = 'pairedEnd fastq file (read1)', required=True)
    parser.add_argument('-o','--outdir', help = 'result directory that all resulting/intermediate files will be stored\n' + \
                                         'will create 1. $resultpath/trimmed\n' + \
                                         '            2. $resultpath/hisat\n'  + \
                                         '            3. $resultpath/bowtie2\n' + \
                                         '            4. $resultpath/mergeBam (all useful result files)\n', required=True)
    parser.add_argument('-x', '--humanIndex', help = 'human bowtie2 index', required=True)
    parser.add_argument('-b','--bedpath', help = 'bed folder for gene counting', required=True)
    parser.add_argument('-s','--splicesite', help = 'splice site file generated by hisat', required=True)
    parser.add_argument('-t','--tRNAindex' , help = 'bowtie2 index for tRNA, for better tRNA counting', required=True)
    parser.add_argument('-r','--rRNAindex' , help = 'bowtie2 index for rRNA, for better rRNA counting', required=True)
    parser.add_argument('-p', '--threads', default=1, type=int, help = 'number of cores to be used for the pipeline (default:1)')
    parser.add_argument('-f', '--adaptors', help = 'fasta file containing adaptors', required=True)
    parser.add_argument('-d', '--strand', choices=['forward','reverse','both'], help = 'strandeness of RNA-seq, can be <forward|reverse|both> default: forward', required=True)
    args = parser.parse_args()
    return args

def strand2int(strandeness):
    if strandeness == 'forward' :
        strand = 0
    elif strandeness == 'reverse':
        strand = 1
    elif strandeness == 'both':
        strand = 2
    return strand

def main():
    programname = sys.argv[0]
    args = getopt()
    fastqFile = args.fastq
    resultpath = args.outdir
    cores = args.threads
    humanIndex = args.humanIndex
    bedpath = args.bedpath
    tRNA_index = args.tRNAindex
    rRNA_index = args.rRNAindex
    adaptors = args.adaptors
    spliceFile = args.splicesite
    strand = strand2int(args.strand)
    sys.stderr.write('Using %s strand for %s\n' %(args.strand,fastqFile) )
    programSequenceControl(fastqFile,resultpath,cores,humanIndex, bedpath,tRNA_index, rRNA_index, adaptors,spliceFile,strand)
    return 0

if __name__ == '__main__':
    main()
