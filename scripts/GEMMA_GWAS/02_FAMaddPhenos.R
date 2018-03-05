#Nicole E Soltis
#02/21/18
#02_FAMaddPhenos

#----------------------------------------------------------------------------
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
Phenos <- read.csv("data/allreadsGWAS/BO5.10/01_prepFiles/lsmeans_zscale_allreads.csv")
Phenos <- Phenos[,-c(1)]

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/B05_GEMMA/01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#col2 = V2 = Isolate
names(Phenos)[1] <- "V2"

Phenos_match <- Phenos[ order(Phenos$V2), ]
myFAM_match <- myFAM
myFAM_match$delete <- c(1:96)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#so 1.01.12 and 1.02.13 need dummy variables
#add an empty variable 1.02.12 to Phenos
for (y in 2:length(Phenos_match)){
  Phenos_match[97,y] <- mean(Phenos_match[,y], na.rm=TRUE)
}
levels(Phenos_match$V2) <- c(levels(Phenos_match$V2), "1.01.12")
Phenos_match[97,1] <- "1.01.12"
#reorder phenos
Phenos_match <- Phenos_match[c(1:5,97,6:11,13:96),]

#now add Phenos_match onto myFAM_match
myFAM_match2 <- cbind(myFAM_match, Phenos_match)
myFAM_match2 <- myFAM_match2[order(myFAM_match$delete),]
myFAM_match2 <- myFAM_match2[,c(1:5,9:length(myFAM_match2))]

Sys.time()
write.table(myFAM_match2, "data/B05_GEMMA/02_GEMMA/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
Sys.time()