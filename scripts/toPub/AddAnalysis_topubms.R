#02142019
#Nicole E Soltis
#analysis to add for publication (post Dissertation)

#---------------------------------------------------------------------------------------------------------------------
#look at permut summary data from At
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
fullpeaks <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
hipeaks <- fullpeaks[fullpeaks$numGenes > 10,]
fullsumm <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/HotspotSumm_5xRand.csv")

i  <- 4
myranddat <- read.csv(paste0("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=''))
#summary info of p-values
quantile(myranddat$p_score, c(0.00, 0.01, 0.05, 0.6))

#          0%           1%           5% 
#3.461793e-08 6.564719e-06 3.014453e-05
#4.923707e-08 5.760712e-06 2.742912e-05 
#4.211615e-08 6.002506e-06 2.886701e-05
#5.873342e-08 5.530269e-06 2.790771e-05 
#2.926053e-08 6.571152e-06 3.071604e-05
mean(c(3.014453e-05, 2.742912e-05, 2.886701e-05, 2.790771e-05, 3.071604e-05))

#real data
mydat <- read.csv("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes.csv")
quantile(mydat$p_score, c(0.00, 0.01, 0.05, 0.6))

#check this for top1 as well
mydat100 <- read.csv("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top100SNPsample_genes.csv")
hi100 <- mydat100[mydat100$p_score < 2.901288e-05,]
hisumm <- as.data.frame(table(hi100$Gene))

#---------------------------------------------------------------------------------------------------------------
#look at permut summary data for Bc
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
hipeaks <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/PeaksOver5.csv")
fullpeaks <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
fullsumm <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/HotspotSumm_5xRand.csv")

i  <- 5
myranddat <- read.csv(paste0("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=''))
#summary info of p-values
quantile(myranddat$p_score, c(0.00, 0.01, 0.05))

#          0%           1%           5% 
#1.530313e-07 5.425669e-06 2.172780e-05 
#2.092033e-07 4.801529e-06 1.864122e-05 
#1.835376e-07 5.526258e-06 1.950385e-05 
#2.949698e-08 4.325383e-06 2.021054e-05 
#1.034973e-07 5.034284e-06 1.795980e-05

mean(c(2.172780e-05, 1.864122e-05, 1.950385e-05, 2.021054e-05, 1.795980e-05))

mydat100 <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
hi100 <- mydat100[mydat100$p_score < 1.960864e-05,]
hisumm <- as.data.frame(table(hi100$pheno))
