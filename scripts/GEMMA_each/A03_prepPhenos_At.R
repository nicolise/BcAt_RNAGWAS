
#Nicole E Soltis 
#03/28/18
#add Phenos to FAM file for GEMMA
#----------------------------------------------------------------------
rm(list=ls())

#pipeline note:
#1. copy files from GEMMA_eachAt_Bc/01_PLINK

setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
IsoNames <- read.csv("data/Vivian_Bc/IsolateKey_Vivian.csv")
MyReads <- read.csv("data/Vivian_At/data from Wei 20160321/S2_ModelLsMeanAdjusted.csv")
MyReads <- MyReads[,-c(1)]

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/GEMMA_eachBc_At/01_PLINK/binMAF20NA10.fam")

#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.
## start here

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")

#col2 = V2 = Isolate
names(Phenos)[2] <- "V2"

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#remove Control, 1.02.13, 94.4
##check these
Phenos_match <- Phenos_match[c(2:12,14:59,61:98),]

#double check now and save original order for FAM
myFAM_match$delete <- c(1:95)
setdiff(myFAM_match$V2,Phenos_match$V2) #check: is zero
setdiff(Phenos_match$V2,myFAM_match$V2) #check: is zero
intersect(myFAM_match$V2,Phenos_match$V2)

#match and return to original order
myFAM_match2 <- merge(myFAM_match, Phenos_match, by="V2")
myFAM_match2 <- myFAM_match2[order(myFAM_match2$delete),]
#remove column "delete" and dummy phenotype (V6) and reorder column V1, V2
myFAM_match2 <- myFAM_match2[,c(2,1,3:5,9:ncol(myFAM_match2))]

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
write.csv(myFAM_match2, "B05_GEMMA/C03_runGEMMA/binMAF20NA10_fam.csv")
write.table(myFAM_match2, "B05_GEMMA/C03_runGEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)
#if randomizing, copy _allphenos.fam to _rand.fam and also copy 01/.bed and 01/.bim to 02/_rand.
