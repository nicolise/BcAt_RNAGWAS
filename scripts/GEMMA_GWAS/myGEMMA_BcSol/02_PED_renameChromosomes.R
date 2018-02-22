#reformat bigRR output data
#Nicole E Soltis

#--------------------------------------------------------
rm(list=ls())
library(tidyr)
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

myMAP <- read.csv("GEMMA_files/02_csvPrep/fulldata/MAP_dpcharMAF20NA10.csv")
myMAPbk <- myMAP

#split chromosome and segment
myDat <- myMAP
names(myDat)

#let's try making the chrom.seg integers so that R isn't confused
unique(myDat$X.CHROM)
myDat$X.CHROM.F <- as.factor(myDat$X.CHROM)
unique(myDat$X.CHROM.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

myDat$X.CHROM.Int <- recode.vars$newvals[match(myDat$X.CHROM.F, recode.vars$OGvals)]
unique(myDat$X.CHROM.Int)

myMAP <- myDat
myMAP <- myMAP[,c("X.CHROM.Int", "POS")]

#add an arbitrary SNPid column
myMAP$SNPid <- paste("SNP", myMAP$X.CHROM.Int, myMAP$POS, sep=".")

myMAP <- myMAP[,c(1, 3, 2)]

library("caroline")
write.delim(myMAP, "GEMMA_files/02_csvPrep/fulldata/MAP_dpcharMAF20NA10_55chrom.tsv", quote = FALSE, col.names = F, row.names = FALSE, sep = "\t")

myPED <- myPED[,-c(1)]
write.delim(myPED, "GEMMA_files/02_csvPrep/fulldata/PED_dpcharMAF20NA10_55chrom.tsv", quote = FALSE, col.names = F, row.names = FALSE, sep = "\t")
