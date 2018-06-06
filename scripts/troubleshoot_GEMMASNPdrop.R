#Nicole E Soltis
#05/18/18

#troubleshoot dropped SNPs B05.10 GEMMA
#----------------------------------------------------------------
#using same genotype input from B05.10 BcAtGWAS bigRR and GEMMA
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#use same SNP set from B05.10 bigRR
mySNPs <- read.csv("allreadsGWAS/BO5.10/01_prepFiles/hp_binMAF20_20NA.csv")

#Chrom 16:18
mySNPs.sub <- mySNPs[mySNPs$Chrom>15,]
#on Chrom 16 only want mySNPs$Pos > 343603

#Index
mySNPs.sub <- mySNPs.sub[with(mySNPs.sub, order(Chrom, Pos)), ]

#Make plotting variables
mySNPs.sub$Index = NA
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mySNPs.sub$Chrom)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    mySNPs.sub[mySNPs.sub$Chrom==i, ]$Index=mySNPs.sub[mySNPs.sub$Chrom==i, ]$Pos
    #for all other Chromomosomes: 
  }	else {
    #lastbase for Chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of Chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(mySNPs.sub,mySNPs.sub$Chrom==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    mySNPs.sub[mySNPs.sub$Chrom==i, ]$Index=mySNPs.sub[mySNPs.sub$Chrom==i, ]$Pos+lastbase
  }
}

#keep only Chr16 past 343,604
mySNPs.sub <- mySNPs.sub[mySNPs.sub$Index > 343604,]

#0, 1, NA
mySNPs.sub$NAcount <- apply(mySNPs.sub[,4:100], 1, function(x) sum(is.na(x)))

#up to 20.83% missing
mySNPs.sub$missing <- mySNPs.sub$NAcount * 100 / 96

#MAF
#this sums how many isolates with allele 1 per SNP
mySNPs.sub$FreqAl1 <- apply(mySNPs.sub[,4:100], 1, function(x) sum(x==1, na.rm=TRUE))
mySNPs.sub$FreqAl0 <- apply(mySNPs.sub[,4:100], 1, function(x) sum(x==0, na.rm=TRUE))
#minor allele count: if allele 1 most common, MAC is 96 - allele 1. 
mySNPs.sub$MAC <- ifelse(mySNPs.sub$FreqAl1 > 50, mySNPs.sub$FreqAl0, mySNPs.sub$FreqAl1)
mySNPs.sub$MajAC <- ifelse(mySNPs.sub$FreqAl1 < 50, mySNPs.sub$FreqAl0, mySNPs.sub$FreqAl1)
mySNPs.sub$TotAlleles <- mySNPs.sub$MAC + mySNPs.sub$MajAC
table(mySNPs.sub$TotAlleles)
hist(mySNPs.sub$TotAlleles)
#is this accurate?
mySNPs.sub$MAF <- mySNPs.sub$MAC / (mySNPs.sub$TotAlleles + mySNPs.sub$missing)
hist(mySNPs.sub$missing)

mySNPs.drop <- mySNPs.sub[mySNPs.sub$missing > 10,]
