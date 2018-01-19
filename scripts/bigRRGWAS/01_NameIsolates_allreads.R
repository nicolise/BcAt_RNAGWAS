#092217
#Nicole E Soltis
#01_NameIsolates.R
#short script to match phenotypes of Isolates with IDs 1-97 with actual Isolate names for matching to SNP data.

#-----------------------------------------------------------------
rm(list=ls())
#load data
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
setwd("~/Projects/BcAt_RNAGWAS/data")
IsoNames <- read.csv("IsolateKey_Vivian.csv")
MyReads <- read.csv("Vivian_Bc/result.lsm.csv")
MyReads <- MyReads[,-c(1)]

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(1,2,9270,3:9269)]

#each plant x transcript is a unique phenotype
#split by host
splitlist <- split(MyReads, MyReads$HostGenotype)
unique(MyReads$HostGenotype)
MyReads.coi1 <- splitlist[[1]]
MyReads.col0 <- splitlist[[2]]
MyReads.npr1 <- splitlist[[3]]

#clean up columns
MyReads.col0 <- MyReads.col0[,c(3,2,4:9270)]
names(MyReads.col0)[1] <- "Isolate"
write.csv(MyReads.col0, "allreadsGWAS/01_prepFiles/col0_allreads.csv")
