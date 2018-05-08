#Nicole E Soltis
#05/04/18 Star Wars day!
#----------------------------------------------------------
rm(list=ls())

setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
Phenos <- read.csv("LSMeanCamLes4Map_FIN.csv")

#phenotype looks fine
hist(Phenos$Col0.Les)
hist(Phenos$coi1.Les)

#try: GEMMA with z-scaled Col0
#GEMMA with no kmat, kmat 1, kmat 2 
#GEMMA with dropping genos/ SNPs from other Phenos
#GEMMA with increased burn in time?

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
myFAM_match2 <- read.csv("B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_fam.csv")

#all phenos NA for 01.01.12, no others
myFAM_col0 <- myFAM_match2[,c(1:8)]
#z-scaled lesion
myFAM_col0$Col0.Les.z <- scale(myFAM_col0$Col0.Les, center = TRUE, scale = TRUE)
hist(myFAM_col0$Col0.Les.z)
write.table(myFAM_col0, "B05_GEMMA_les/troubleshoot_Col0/D_02_GEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)

