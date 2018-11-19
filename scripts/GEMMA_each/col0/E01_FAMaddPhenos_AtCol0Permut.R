
#Nicole E Soltis 
#03/28/18
#add Phenos to FAM file for GEMMA
#----------------------------------------------------------------------
#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
IsoNames <- read.csv("data/Vivian_Bc/IsolateKey_Vivian.csv")
MyReads <- read.csv("data/Vivian_At/data from Wei 20160321/S2_ModelLsMeanAdjusted.csv")
MyReads <- MyReads[,-c(1)]
#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(23960,2,4:23959)]
names(MyReads)[1] <- "Isolate"

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
#*.bed, *.bim, and *.fam are copied from B05_GEMMA_les/E_01_PLINK -- same as in B05_GEMMA_Bc but with 01.01.12 removed so fewer SNPs will be dropped from GEMMA
myFAM <- read.table("data/GEMMA_eachBc_At/01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#first, split MyReads by plant accession. Then generate matching FAM files for each.
myread_col0 <- MyReads[MyReads$HostGenotype=="col.0",]
Phenos <- myread_col0
#col2 = V2 = Isolate
Phenos <- Phenos[,c(1,3:length(Phenos))]
names(Phenos)[1] <- "V2"

#randomize each Phenotype!
Phenos_rand <- Phenos
for (i in 2:23957){
  Phenos_rand[,i] <- sample(Phenos_rand[,i])
}
#select columns
Phenos_match <- Phenos_rand[ order(Phenos_rand$V2), ]
#remove non-genotyped 01.02.13 from Phenos
Phenos_match <- Phenos_match[!Phenos_match$V2 =="1.02.13",]
myFAM_match <- myFAM
myFAM_match$delete <- c(1:95)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]

## check that these are 0, 0, 95
setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#now add Phenos_match onto myFAM_match
myFAM_match2 <- merge(myFAM_match, Phenos_match, by="V2")
myFAM_match2 <- myFAM_match2[order(myFAM_match2$delete),]
#remove dummy phenotype (column 6)
#and reorder V1:V5
myFAM_match2 <- myFAM_match2[,c(2,1,3:5,8:length(myFAM_match2))]

randFAM5 <- myFAM_match2
##be sure to move new *.fam into the correct directory!
#could run this to Data but the drive is pretty full
#trying to external
#setwd("/media/nesoltis/Data/Kliebenstein/Soltis/BcAt_RNAGWAS")
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS")
Sys.time()
write.table(randFAM1, "GEMMA_eachBc_At/col0/02_GEMMA_RAND01/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
write.table(randFAM2, "GEMMA_eachBc_At/col0/02_GEMMA_RAND02/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
write.table(randFAM3, "GEMMA_eachBc_At/col0/02_GEMMA_RAND03/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
write.table(randFAM4, "GEMMA_eachBc_At/col0/02_GEMMA_RAND04/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
write.table(randFAM5, "GEMMA_eachBc_At/col0/02_GEMMA_RAND05/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
Sys.time()
#----------------------------------------------
#location of original (actual) data:
##be sure to move new *.fam into the correct directory!
Sys.time()
write.table(myFAM_match2, "data/GEMMA_eachBc_At/02_GEMMA/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
myFAM_check <- read.table("data/GEMMA_eachBc_At/02_GEMMA/binMAF20NA10.fam")
Sys.time()