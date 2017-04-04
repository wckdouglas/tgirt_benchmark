args<-commandArgs(TRUE)
file_read=args[1]
file_write=args[2]
rRNA<-c("gi|23898|emb|X12811.1|", "gi|555853|gb|U13369.1|HSU13369")
chr<-c(1:22,"X","Y","MT")
conR=file(file_read,open="r")
conW=file(file_write,open="a")
reads1<-readLines(conR,n=1,warn=FALSE)
l1<-strsplit(reads1,'\t')
reads2<-readLines(conR,n=1,warn=FALSE)
l2<-strsplit(reads2,'\t')
# pre-set the FLAG, to indicate if reads mapped to rRNA
if (l1[[1]][3] %in% rRNA) {
  FLAG1=TRUE
}else{
  FLAG1=FALSE
  if (l1[[1]][3] %in% chr){
    chr1=TRUE
  }else{
    chr1=FALSE
  }
}
TLEN1<-abs(as.numeric(l1[[1]][9]))
ID1<-l1[[1]][1]
repeat {
  reads3<-readLines(conR,n=1,warn=FALSE)
  # check if read to the end of the file, if so, write the last two previous reads and stop
  if (length(reads3)==0)
  {
    writeLines(reads1,conW) 
    writeLines(reads2,conW) 
    break
  }
  # if not end, read new set of reads, and assign TLEN, ID
  l3<-strsplit(reads3,'\t')
  reads4<-readLines(conR,n=1,warn=FALSE)
  l4<-strsplit(reads4,'\t')
  if (l3[[1]][3] %in% rRNA) {
    FLAG2=TRUE
  }else{
    FLAG2=FALSE
    if (l3[[1]][3] %in% chr){
    chr2=TRUE
    }else{
    chr2=FALSE
  }
  }
  TLEN2<-abs(as.numeric(l3[[1]][9]))
  ID2<-l3[[1]][1]
  if (ID2!=ID1){
    # if the IDs are different, write the pair reads 1 to file, 
    # and save read pair 2 to read pair1, reset the FLAG
    writeLines(reads1,conW) 
    writeLines(reads2,conW) 
    reads1<-reads3
    reads2<-reads4
    TLEN1<-TLEN2
    ID1<-ID2
    FLAG1<-FLAG2
    chr1<-chr2
  }else{
    if (TLEN2<TLEN1){
      # if the new set of reads has shorter TLEN
      # replace set1 with set2, reset FLAG
      reads1<-reads3
      reads2<-reads4
      l1<-l3
      l2<-l4
      TLEN1<-TLEN2
      ID1<-ID2
      FLAG1<-FLAG2
      chr1<-chr2
    }else if (TLEN2==TLEN1){
      # if the two sets of reads have the same TLEN
      # check FLAG, keep the rRNA reads, discard the other
      # if all are rRNA, keep set1
      if (!FLAG1 & (FLAG2 |(!chr1 & chr2))) {
        reads1<-reads3
        reads2<-reads4
        l1<-l3
        l2<-l4
        TLEN1<-TLEN2
        ID1<-ID2
        FLAG1<-FLAG2
        chr1<-chr2
      }
    }
  }
}
close(conR)
close(conW)
