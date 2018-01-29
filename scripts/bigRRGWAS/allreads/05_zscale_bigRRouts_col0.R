#Nicole E Soltis
#012218

#z-scale bigRR outputs from each transcript to directly compare across reads
#------------------------------------------------------------
library(tidyr)

#read in bigRR output file
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcAt_RNAGWAS/data/")
#Import data
#first stick the two output files together
mydat01 <- read.csv("allreadsGWAS/03_bigRRouts/col0_bigRR_MAF20_011818.csv")
mydat02 <- read.csv("allreadsGWAS/03_bigRRouts/col0_bigRR_MAF20_011918.csv")
#remove empty columns
mydat01 <- mydat01[,-c(1)]
mydat02 <- mydat02[,-c(1)]
head(mydat01)[1:5]
head(mydat02)[1:5]

#so mydat01 only ran for the first 2.2 chromosomes
#work with just mydat02 for now, only got 226 of the reads.

#split chromosome and segment
mydat <- separate (mydat02, outpt.HEM, into = c("Chrom", "Segment", "Pos") )
#double check
unique(mydat$Chrom)
unique(mydat$Segment)

#Reformat Chromosomes and Positions
mydat$Chrom <- as.numeric(as.character(mydat$Chrom))
mydat$Segment <- as.numeric(as.character(mydat$Segment))
mydat$Pos <- as.numeric(as.character(mydat$Pos))

#sort dataframe rows in order of Chrom, then Pos
mydat <- mydat[with(mydat, order(Chrom, Segment, Pos)), ]

#now z-scale the phenotypes
mydat_z <- mydat[,1:3]
#4:228
for (i in c(4:228)){
  mydat_z[,i] <- scale(mydat[,i], center = TRUE, scale = TRUE)
  names(mydat_z)[i] <- print(names(mydat)[i])
}

#now build new dataframe for synteny thing
#add variable chrom.seg.pos
#columns: Chromosome, Segment, Pos, chrom.seg.pos, Transcript
mydat_z$chrom.seg.pos <- paste(mydat_z$Chrom, mydat_z$Segment, mydat_z$Pos, sep=".")
mydat_z <- mydat_z[,c(1:3, 229, 4:228)]

#keep only top SNP per transcript for cis/trans analysis
mydat_cistrans <- as.data.frame(names(mydat_z)[5:10])
names(mydat_cistrans)[1] <- "Transcript"

#start the dataframe
#first phenotype is 5
df <- mydat_z[mydat_z[,5]==max(abs(mydat_z[,5])),]
df$Transcript <- names(df)[5]
df$maxFX <- df[,5]
df <- df[,c(1:4,230, 231)]

for (i in c(6:229)){
#extract the row
blah <- mydat_z[abs(mydat_z[,i])==max(abs(mydat_z[,i])),]
blah$Transcript <- names(blah)[i]
blah$maxFX <- blah[,i]
blah <- blah[,c(1:4,230, 231)]
df <- rbind(df, blah)
}

#now extract positional information per SNP
myreadPos <- as.data.frame(names(mydat_z)[4:228])
names(myreadPos)[1] <- "readname"

#very bad lazy regex here
myreadPos <- df
names(myreadPos)[1] <- "SNP.Chrom"
names(myreadPos)[2] <- "SNP.Seg"
names(myreadPos)[3] <- "SNP.Pos"
names(myreadPos)[4] <- "SNP.c.s.p"
myreadPos$read.Chrom <- substring(myreadPos$Transcript, regexpr("Bcin", myreadPos$Transcript) +4)
myreadPos$read.Chrom <- as.numeric(sapply(strsplit(myreadPos$read.Chrom, "g"), "[", 1))

myreadPos$read.Pos <- substring(myreadPos$Transcript, regexpr("g", myreadPos$Transcript) + 1)
myreadPos$read.Pos <- sapply(strsplit(myreadPos$read.Pos, ".1.HEM"), "[", 1)
myreadPos$read.Pos <- sapply(strsplit(myreadPos$read.Pos, ".2.HEM"), "[", 1)
myreadPos$read.Pos <- as.numeric(myreadPos$read.Pos)

#include read 1 or 2 for a given SNP
myreadPos$read.num <- substring(myreadPos$Transcript, regexpr("g", myreadPos$Transcript) + 1)
myreadPos$read.num <- as.numeric(as.character(substring(myreadPos$read.num, regexpr("g", myreadPos$read.num), regexpr(".HEM", myreadPos$read.num) -1)))

#have no contig info for reads so... for now ignoring contig for SNP and just going Chrom/ Pos
#may need to go back to old SlBcGWAS bigRRout script to grab Chr.Contig indexing

#must do this twice: once for SNP.Index and once for read.Index
#first order data by SNP.Chrom
myreadPos <- myreadPos[order(myreadPos$SNP.Chrom),]
#Make plotting variables
myreadPos$SNP.Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myreadPos$SNP.Chrom)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Index=myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myreadPos,myreadPos$SNP.Chrom==i-1)$SNP.Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Index=myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Pos+lastbase
  }
  ticks=c(ticks, myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Index[floor(length(myreadPos[myreadPos$SNP.Chrom==i, ]$SNP.Index)/2)+1])
}

#now order data by read.Chrom
myreadPos <- myreadPos[order(myreadPos$read.Chrom),]
#Make plotting variables
myreadPos$read.Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myreadPos$read.Chrom)) {
  print(i)
  #for chromosome 1
  if (i==11) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myreadPos[myreadPos$read.Chrom==i, ]$read.Index=myreadPos[myreadPos$read.Chrom==i, ]$read.Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myreadPos,myreadPos$read.Chrom==i-1)$read.Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myreadPos[myreadPos$read.Chrom==i, ]$read.Index=myreadPos[myreadPos$read.Chrom==i, ]$read.Pos+lastbase
  }
}

#plot it
library(ggplot2)
ggplot(myreadPos, aes(x=SNP.Index, y=read.Index)) + geom_point(aes(color=read.Chrom, size=SNP.Chrom)) 
  

write.csv(HEM.plotdata, "BcBOT_MAF20_HEM.PlotFormat.csv") 
write.csv(HEMthresh, "BcBOT_MAF20_HEM.Thresh.csv")
#read in to 06_plots