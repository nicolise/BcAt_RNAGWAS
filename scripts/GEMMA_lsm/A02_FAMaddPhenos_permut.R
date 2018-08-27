#Nicole E Soltis
#02/21/18
#02_FAMaddPhenos

#----------------------------------------------------------------------------
rm(list=ls())

#GEMMA run script is /home/nesoltis/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA/runGEMMA_allphenos.sh 

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
Phenos <- read.csv("data/allreads_bigRR/B05.10/01_prepFiles/lsmeans_zscale_allreads.csv")
Phenos <- Phenos[,-c(1)]

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam

#*.bed, *.bim, and *.fam are copied from B05_GEMMA_les/E_01_PLINK -- same as in B05_GEMMA_Bc but with 01.01.12 removed so fewer SNPs will be dropped from GEMMA
myFAM <- read.table("data/B05_GEMMA_Bc/B_permut/01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#col2 = V2 = Isolate
names(Phenos)[1] <- "V2"

#remove Phenotype 01.02.13 because no genotype info
Phenos_match <- Phenos_match[-c(12),]

#randomize each Phenotype!
Phenos_rand <- Phenos
for (i in 2:9268){
Phenos_rand[,i] <- sample(Phenos_rand[,i])
}
#select columns
Phenos_match <- Phenos_rand[ order(Phenos_rand$V2), ]
myFAM_match <- myFAM
myFAM_match$delete <- c(1:95)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#now add Phenos_match onto myFAM_match
myFAM_match2 <- merge(myFAM_match, Phenos_match, by = "V2")
myFAM_match2 <- myFAM_match2[order(myFAM_match$delete),]

#correctly order V1:V5, remove dummy phenotype V6
myFAM_match2 <- myFAM_match2[,c(2,1,3:5,8:length(myFAM_match2))]

Sys.time()
write.table(myFAM_match2, "data/B05_GEMMA_Bc/B_permut/02_GEMMA/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
myFAM_check <- read.table("data/B05_GEMMA_Bc/B_permut/02_GEMMA/binMAF20NA10.fam")
Sys.time()