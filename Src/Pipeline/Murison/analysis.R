# Set this to the current working directory
setwd("/media/mentat/HiC/LowC/HinDIII/100K/aligned/")

library(GenomicRanges)

DPNFile<-"/media/mentat/HiC/hg19_HindIII_new.txt"
con=file(DPNFile,open="r")
line=readLines(con) 
long=length(line)

i<-1
splitcurrent<-strsplit(line[i], "\\ ")
chromosome<-splitcurrent[[1]][1]
mymatrix<-matrix(ncol=3, nrow=length(splitcurrent[[1]])-1)
mymatrix[,1]<-rep(chromosome, nrow(mymatrix))
for(j in 1:nrow(mymatrix)) {mymatrix[j,2]<-splitcurrent[[1]][j+1]}
mymatrix[,3]<-as.numeric(mymatrix[,2])+3

HINSites<-GRanges(seqnames=as.character(mymatrix[,1]), IRanges(start=as.numeric(mymatrix[,2]), end=as.numeric(mymatrix[,3])))

for(i in 2:length(line)) {
  
  
  splitcurrent<-strsplit(line[i], "\\ ")
  chromosome<-splitcurrent[[1]][1]
  print(chromosome)
  mymatrix<-matrix(ncol=3, nrow=length(splitcurrent[[1]])-1)
  mymatrix[,1]<-rep(chromosome, nrow(mymatrix))
  for(j in 1:nrow(mymatrix)) {mymatrix[j,2]<-splitcurrent[[1]][j+1]}
  mymatrix[,3]<-as.numeric(mymatrix[,2])+3
  
  HINSites<-c(HINSites, GRanges(seqnames=as.character(mymatrix[,1]), IRanges(start=as.numeric(mymatrix[,2]), end=as.numeric(mymatrix[,3]))))
}

setwd("/media/mentat/Prostate/LowC/PCa33173/")

reads<-read.table("merged_nodups.txt", sep=" ", stringsAsFactors = FALSE)
colnames(reads)<-c("str1","chr1","pos1","frag1","str2","chr2","pos2","frag2","mapq1","cigar1","sequence1","mapq2","cigar2","sequence2","readname1","readname2")

### Minimum MapQ of 3
filterme<-unique(c(which(reads$mapq1<3), which(reads$mapq2<3)))
reads_q3<-reads[-filterme,]


### Remove those reads mapping to the same fragment

remfilt<-which((reads_q3$frag1-reads_q3$frag2)==0)
length(remfilt)
reads_q3_DiffFrag<-reads_q3[-remfilt,]

### Remove those reads > 5kb away from a fragment

F1<-GRanges(seqnames=as.character(reads_q3_DiffFrag$chr1), IRanges(start=as.numeric(reads_q3_DiffFrag$pos1), end=as.numeric(reads_q3_DiffFrag$pos1+nchar(reads_q3_DiffFrag$sequence1)-1)))
F2<-GRanges(seqnames=as.character(reads_q3_DiffFrag$chr2), IRanges(start=as.numeric(reads_q3_DiffFrag$pos2), end=as.numeric(reads_q3_DiffFrag$pos2+nchar(reads_q3_DiffFrag$sequence2)-1)))

nearestF1<-distanceToNearest(F1, HINSites)
nearestF2<-distanceToNearest(F2, HINSites)

Dm<-cbind(mcols(nearestF1)[,1], mcols(nearestF2)[,1])

mindist<-apply(Dm, MARGIN = 1, min)
reads_q3_DiffFrag_Within5Kb<-reads_q3_DiffFrag[which(mindist<5000),]

### I have yet to implement the filtering based upon sequencing bias

### Output Filtered Read

write.table(reads_q3_DiffFrag_Within5Kb[,-ncol(reads_q3_DiffFrag_Within5Kb)], "merged_nodups_q3_DiffFrag_Within5k.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)



### Minimum MapQ of 1
filterme<-unique(c(which(reads$mapq1<1), which(reads$mapq2<1)))
reads_q1<-reads#[-filterme,]


### Remove those reads mapping to the same fragment

remfilt<-which((reads_q1$frag1-reads_q1$frag2)==0)
length(remfilt)
reads_q1_DiffFrag<-reads_q1[-remfilt,]

### Remove those reads > 5kb away from a fragment

F1<-GRanges(seqnames=as.character(reads_q1_DiffFrag$chr1), IRanges(start=as.numeric(reads_q1_DiffFrag$pos1), end=as.numeric(reads_q1_DiffFrag$pos1+nchar(reads_q1_DiffFrag$sequence1)-1)))
F2<-GRanges(seqnames=as.character(reads_q1_DiffFrag$chr2), IRanges(start=as.numeric(reads_q1_DiffFrag$pos2), end=as.numeric(reads_q1_DiffFrag$pos2+nchar(reads_q1_DiffFrag$sequence2)-1)))

nearestF1<-distanceToNearest(F1, HINSites)
nearestF2<-distanceToNearest(F2, HINSites)

Dm<-cbind(mcols(nearestF1)[,1], mcols(nearestF2)[,1])

mindist<-apply(Dm, MARGIN = 1, min)
reads_q1_DiffFrag_Within5Kb<-reads_q1_DiffFrag[which(mindist<5000),]

### I have yet to implement the filtering based upon sequencing bias

### Output Filtered Read

write.table(reads_q1_DiffFrag_Within5Kb[,-ncol(reads_q1_DiffFrag_Within5Kb)], "merged_nodups_q1_DiffFrag_Within5k.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

### Take exported file and use Use juicer pre to turn it into a .hic file ( use command "java -jar juicerXXX.jar pre merged_nodups_q1_DiffFrag_Within5k.txt merged_nodups_q1_DiffFrag_Within5k.hic hg19")

### From exported file use juicer dump to output per-chromosome normalized interaction frequencies at the desired resolution e.g.  java -jar juicer_tools.1.8.9_jcuda.0.8.jar dump observed KR merged_nodups_q1_DiffFrag_Within5k.hic 1 1 BP 100000 chr1_Filtered_KR_100kb.txt
# Explained bwlow
# java -jar juicer_tools.1.8.9_jcuda.0.8.jar dump observed KR PATOTOINPUT/merged_nodups_q1_DiffFrag_Within5k.hic \ 
# 1 1 \ #Chromsome A and B to compute interactions over
# BP 100000 \ # BP orientation, and 100kb windows
# PATHTOOUTPUT/chr1_Filtered_KR_100kb.txt

### Juicer outputs these as small files only listing regions where any interactions are observed. For calling TADs I have used the Insulation Score (IS) Method. 
# I either need to add to this a filter whereby I remove centromeric regions prior to running IS or filter them out afterwards... 

# Now I put these into a format which can be read by Insulation_Score. I ended up rerunning and changing these parameters to 500kb resolution

options(scipen=999)
resolution<-100000

for(j in c(1:22, "X", "Y")) {
  
  ### Generate the filename and read
  chrname<-j
  chromosome<-paste0("chr", chrname,"_Filtered_KR_100kb.txt") # Adjust as needed
  KR<-read.table(chromosome, sep="\t", stringsAsFactors = FALSE)
  
  ### Get the highest bin number directly from the matrix and use resolution to create windows
  maxnum<-max(c(KR[,2],KR[,1]))
  matcols<-seq(0,maxnum, resolution)
  
  ### Insulation Score also has some annoyingly specific requirements for how the column and rownames of the interaction matrix are set up
  ### They have to follow the format binXXX|hgQQ|chrAA:YYYYY-ZZZZZ
  ### Here we generate the column names and make the contact matrix
  enw1<-paste0("bin", 1:length(matcols))
  #enw1<-enw1[-length(enw1)]
  enw2<-paste0(enw1, "|hg19|chr",chrname)
  matcolsplusone<-c(matcols, matcols[length(matcols)]+resolution)
  for(i in 1:length(enw2)) {
    enw2[i]<-paste0(enw2[i], ":", matcols[i], "-", matcolsplusone[i+1])
  }
  ContactMatrix<-matrix(ncol=length(enw2), nrow=length(enw2))
  rownames(ContactMatrix)<-matcols
  colnames(ContactMatrix)<-matcols
  
  ### Now we can finally read in the KR balanced terms and slot them into the appropriate bin within the contact matrix
  for(i in 1:nrow(KR)) {
    ContactMatrix[as.character(KR[i,1]),as.character(KR[i,2])]<-KR[i,3]
  }
  
  ### And finally we can format the contact matrix properly and output it
  ContactMatrix[which(is.na(ContactMatrix))]<-0
  rownames(ContactMatrix)<-enw2
  colnames(ContactMatrix)<-enw2
  write.table(ContactMatrix, paste0("chr", chrname,"_Filtered_KR_100kb_ContactMatrix.tsv"), sep="\t", quote=FALSE)
  
}

system("sed -i '1s/^/\t/' chr*_Filtered_KR_100kb_ContactMatrix.tsv")

### Run IS (Insulation Score)
## https://github.com/dekkerlab/crane-nature-2015
## I used 500kb resolution and I set bmoe to 0 
## Comparing to hg19 LNCaP TADs