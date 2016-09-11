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
    command='trimmomatic PE '' +\
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
    split_command = 'python splitBam.py --outprefix=%s --program=hisat2 --inBam=%s' %(outFile,outFile)
    runProcess((split_command,sampleName))
    return outFile

def getID(bamFile,IDpath,sampleName):
    """
        converting unmapped bamfile to ID file
        input:
            bam file
            resultpath
            sample name
    """
    idFile = '%s/%s.id.dat' %(IDpath,sampleName)
    command = "%s/samtools view -F 4 -@ 12 %s " %(samtoolsPath,bamFile)+\
            "| awk {'print $1'} > %s" %idFile
    runProcess((command,sampleName))
    return idFile

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
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    resultfile = '%s/%s' %(resultpath,sampleName)
    command = '%s/bowtie2 --local  --threads %i  ' %(bowtiePath, cores)+ \
            '--time -L 8 -x %s -1 %s -2 %s ' %(index,file1,file2) + \
            '| awk \'$1~"@" || $2==163 || $2==83 ||  $2==99 || $2==147\'' +\
            '| %s/samtools view -@ %i -b - ' %(samtoolsPath,cores) + \
            '| %s/filterSoftClipped -s 0.3 -b 0.4 - ' %(filterSamToolsPath)+ \
            '| tee %s.bam' %(resultfile) +\
            '| %s/python %s/splitBam.py --outprefix=%s --inBam=- --program=bowtie2' %(pythonPath,scriptDir,resultfile)
    runProcess((command,sampleName))
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
    command = '%s/samtools cat  %s %s ' %(samtoolsPath,bam1,bam2)+\
	    '| %s/samtools fixmate -O bam -r - - ' %(samtoolsPath) +\
            '> ' + outBam
    runProcess((command,sampleName))
    return outBam

def bamToBed(bamFile,resultpath,sampleName):
    """
        converting bamfile to bedfile and count the features
        input:
            bam file
            bed file
            result path for storing count files
            sample name
    """
    bedFile = '%s/%s.bed' %(resultpath,sampleName)
    command = "%s/bamToBed -bedpe -mate1 -i %s " %(bedtoolsPath, bamFile) +\
            '| %s/bedpeTobed -i - -m 1000 > %s' %(bedFileToolsPath, bedFile)
    runProcess((command,sampleName))
    return bedFile

def bed2Count(bedfile, geneBed, tRNAbed,resultpath,sampleName,strand):
    tmprDir = resultpath+'/'+sampleName
    makeFolder(tmprDir)
    if strand == 0:
        type = ' -s '
    elif strand == 1:
        type = ' -S '
    elif strand == 2:
        type = ' '
    baseCountFile = "%s/%s.baseCounts" %(resultpath,sampleName)
    command = '%s/bedtools intersect -a %s -b %s -f 0.5 -wb -v ' %(bedtoolsPath,bedfile,tRNAbed) +\
            '| %s/bedtools intersect -a - -b %s -f 0.5 -wb %s ' %(bedtoolsPath, geneBed, type)+\
            "| awk '{print $NF}' " + \
            '| sort --temporary-directory=%s ' %(tmprDir)+\
            '| uniq -c ' +\
            "| awk '{print $2,$1}' OFS='\\t' "\
            ">  %s" %(baseCountFile)
    runProcess((command, sampleName))
    os.rmdir(tmprDir)
    return baseCountFile

def multiBed2Count(multiBed, geneBed, resultpath, sampleName, strand):
    if strand == 0:
        type = ' -s '
    elif strand == 1:
        type = ' -S '
    elif strand == 2:
        type = ' '
    command = '%s/bedtools intersect -a %s -b %s -f 0.8 -wb -bed %s' %(bedtoolsPath,multiBed, geneBed, type) + \
            "| awk '{print $4,$NF}' OFS='\\t' " +\
            "> %s/%s.multiCounts" %(resultpath, sampleName)
    runProcess((command,sampleName))
    return 0

def getUnmappedID(uniqueBam, multiBam, idpath, sampleName, tRNAbed):
    idFile = '%s/%s.id.dat' %(idpath,sampleName)
    command = "%s/samtools cat %s %s " %(samtoolsPath, uniqueBam, multiBam)+ \
            "| %s/samtools view -bF 4 " %(samtoolsPath)+\
            '| %s/bedtools intersect -abam - -b %s -f 0.1 -v ' %(bedtoolsPath, tRNAbed) + \
            "| %s/samtools view " %(samtoolsPath)+\
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
    commands = ["%s/filterFastq -q %s/%s_%dP.fastq.gz -i %s > %s/%s_%dP.fastq" \
            %(fastqToolsPath,fastqPath,sampleName,end,\
            idFile,resultpath,sampleName,end) for end in [1,2]]
    Pool(cores).map(runProcess,[(command, sampleName) for command in commands])
    return 0

def humantRNA(datapath,cores,sampleName,resultpath,tRNA_index,strand):
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
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    if strand == 0:
        type = ' --norc '
    elif strand == 1:
        type = ' --nofw '
    elif strand == 2:
        type = ' '
    idFile = '%s/%s.id.dat' %(resultpath,sampleName)
    command = '%s/bowtie2 --threads %i -L 18  ' %(bowtiePath,cores)+ \
            type + \
            '--gbar 10 -x %s -1 %s -2 %s ' %(tRNA_index, file1,file2) + \
            "| grep -v \'^@\' " + \
            '| awk \'$3!="\*" {print $1}\' ' + \
            '> %s ' %idFile
    runProcess((command, sampleName))
    return idFile

def mappingTRNA(datapath,cores,sampleName,resultpath,tRNA_index,strand):
    """
        Reassign tRNA reads to tRNA species
        and count the species
        input:
            fastq file
            cores
            sample name
            result directory
    """
    start = time.time()
    file1 = datapath + '/' + sampleName + '_1P.fastq'
    file2 = datapath + '/' + sampleName + '_2P.fastq'
    tmprDir = resultpath+'/'+sampleName
    makeFolder(tmprDir)
    if strand == 0:
        type = ' --norc '
    elif strand == 1:
        type = ' --nofw '
    elif strand == 2:
        type = ' '
    command = "%s/bowtie2 --local --threads %i -L 18  " %(bowtiePath,cores) + \
            type + \
            "--gbar 10 -x %s -1 %s -2 %s " %(tRNA_index, file1,file2) +\
            "| awk '$1~\"@\" || $2==163 || $2==83 || $2==99 || $2==147' " +\
            '| %s/samtools view -h@ %i -bq 1 - ' %(samtoolsPath,cores) +\
            '| tee %s/%s.bam ' %(resultpath,sampleName) +\
            '| %s/samtools view -@ %i - ' %(samtoolsPath,cores) +\
            '| cut -f1,3' +\
            "| sort --temporary-directory=%s " %(tmprDir) +\
            '| uniq -c '+\
            "| awk '$1 == 2 {print $3}' " + \
            "| sort --temporary-directory=%s " %(tmprDir) +\
            "| uniq -c " + \
            "| awk '{print $2,$1}' OFS='\\t' " +\
            "| sort --temporary-directory=%s -k2nr " %(tmprDir) + \
            "> %s/%s.tRNA.counts" %(resultpath,sampleName)
    runProcess((command, sampleName))
    os.rmdir(tmprDir)
    return 0

def programSequenceControl(fastqFile,resultpath,cores,humanIndex,genesBed,tRNA_index,tRNAbed,adaptors,spliceFile,strand):
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
    humanTRNAFastqPath = mergedPath + '/humanTRNA'

    # make result directories
    folders = [resultpath,trimResultPath,hisatResultPath,hisatUnmappedPath,hisatMappedPath,
            bowtieResultPath,bowtieMappedPath,mergedPath,mergedBamPath,uniqueBamPath, multiBamPath,
            countPath,humanTRNAFastqPath,alltRNAfastqPath, mergedBedPath, uniqueBedPath, multiBedPath]
    map(makeFolder,folders)

    # trimmomatic
    trimmomatic(fastqFile, trimResultPath,sampleName,cores,adaptors)

    # hisat map
    hisatOutprefix = hisat_pariedEnd( trimResultPath,
            hisatMappedPath,
            sampleName, cores ,
            humanIndex, spliceFile)

    # extract mapped ID
    idFile = getID(hisatOutprefix + '.bam', hisatMappedPath, sampleName)

    # extract unmapped sequence
    fastqRemoveID(sampleName,trimResultPath,hisatUnmappedPath, idFile, cores)

    # local mapped
    bowtieOutprefix = bowtie_pairedEnd(hisatUnmappedPath, bowtieResultPath,
                                        humanIndex, sampleName,cores)

    #merge mapped bams
    multiBam, uniqueBam = [mergeBams(bowtieOutprefix + '.%s.bam' %maptype,
                                        hisatOutprefix + '.%s.bam' %maptype,
                                        path, sampleName,cores)  \
            for maptype, path in zip(['multi','unique'],[multiBamPath,uniqueBamPath])]

    #converting bam file to bed file
    uniqueBed, multiBed = [bamToBed(bam, bedpath, sampleName) \
            for bam, bedpath in zip([uniqueBam,multiBam],[uniqueBedPath,multiBedPath])]

    #count bed
    baseCountFile = bed2Count(uniqueBed, genesBed, tRNAbed, countPath, sampleName, strand)
    multiCount = multiBed2Count(multiBed, genesBed, countPath, sampleName, strand)

    #all unique Reads ID
    idFile = getUnmappedID(uniqueBam, multiBam, alltRNAfastqPath, sampleName, tRNAbed)

    #extract all sequence for mapping tRNA
    fastqRemoveID(sampleName, trimResultPath, alltRNAfastqPath, idFile, cores)

    #mapping tRNA
    idFile = humantRNA(alltRNAfastqPath,cores,sampleName,humanTRNAFastqPath,tRNA_index,strand)

    #extract all sequence for mapping tRNA
    id2Fastq(trimResultPath,sampleName, idFile, humanTRNAFastqPath, cores)

    #final map tRNA
    mappingTRNA(humanTRNAFastqPath,cores,sampleName,countPath,tRNA_index, strand)

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
    parser.add_argument('-b','--genesBed', help = 'human bed file for gene counting', required=True)
    parser.add_argument('-s','--splicesite', help = 'splice site file generated by hisat', required=True)
    parser.add_argument('-r','--tRNAindex' , help = 'bowtie2 index for tRNA, for better tRNA counting', required=True)
    parser.add_argument('-a','--tRNAbed', help = 'tRNA bed file for removing tRNA reads from initial mapping', required=True)
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
    genesBed = args.genesBed
    tRNA_index = args.tRNAindex
    tRNAbed = args.tRNAbed
    adaptors = args.adaptors
    spliceFile = args.splicesite
    strand = strand2int(args.strand)
    sys.stderr.write('Using %s strand for %s\n' %(args.strand,fastqFile) )
    programSequenceControl(fastqFile,resultpath,cores,humanIndex,genesBed,tRNA_index,tRNAbed,adaptors,spliceFile,strand)
    return 0

if __name__ == '__main__':
    main()
