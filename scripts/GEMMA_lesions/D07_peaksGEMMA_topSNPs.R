#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcSolGWAS/BcAt_RNAGWAS")
#first round: just one file at a time. Then convert to loops

#get thresholds here 
mythrs <- read.csv("data/B05_GEMMA_les/D_07_randSUMM/GEMMA_1krand_Thresholds.csv")
mythrs

#columns with pscore: 13 Col0, 16 coi1, 19 npr1, 22 pad3
myGEMMA <- read.csv("data/B05_GEMMA_les/D_06_results/LesionPhenos_allSNPs_MAF20NA10_GEMMA_kmat1_Indexed.csv")
fulldat <- myGEMMA

#get thresholds here 
mythr999 <- mythrs[mythrs$SNPnum == 236,]
mythr99 <- mythrs[mythrs$SNPnum == 2357,]
mythr9999 <- mythrs[mythrs$SNPnum == 24,]

## check which thr
curthrs <- mythr999
thrCol0 <- curthrs[curthrs$pheno==1,6]
thrcoi1 <- curthrs[curthrs$pheno==2,6]
thrnpr1 <- curthrs[curthrs$pheno==3,6]
thrpad3 <- curthrs[curthrs$pheno==4,6]

#select just top SNPs for comparison to bigRR T4
#conditionally replace nonsig values with zero
myGEMMA$X1_Col0.Les_beta[myGEMMA$X1_Col0.Les_pscore > thrCol0] <- 0
myGEMMA$X2_coi1.Les_beta[myGEMMA$X2_coi1.Les_pscore > thrcoi1] <- 0
myGEMMA$X3_npr1.Les_beta[myGEMMA$X3_npr1.Les_pscore > thrnpr1] <- 0
myGEMMA$X4_pad3.Les_beta[myGEMMA$X4_pad3.Les_pscore > thrpad3] <- 0

#remove rows if all 3 = 0 
myGEMMA_2 <- myGEMMA[!(myGEMMA$X1_Col0.Les_beta==0 & myGEMMA$X2_coi1.Les_beta==0 & myGEMMA$X3_npr1.Les_beta==0 & myGEMMA$X4_pad3.Les_beta==0),]

myGEMMA.col0 <- myGEMMA[,c(3:11,30,12:14)]
myGEMMA.col0 <- myGEMMA.col0[!myGEMMA.col0$X1_Col0.Les_beta==0,]
myGEMMA.coi1 <- myGEMMA[,c(3:11,30,15:17)]
myGEMMA.coi1 <- myGEMMA.coi1[!myGEMMA.coi1$X2_coi1.Les_beta==0,]
myGEMMA.npr1 <- myGEMMA[,c(3:11,30,18:20)]
myGEMMA.npr1 <- myGEMMA.npr1[!myGEMMA.npr1$X3_npr1.Les_beta==0,]
myGEMMA.pad3 <- myGEMMA[,c(3:11,30,21:23)]
myGEMMA.pad3 <- myGEMMA.pad3[!myGEMMA.pad3$X4_pad3.Les_beta==0,]

##check which files
write.csv(myGEMMA.pad3, "data/B05_GEMMA_les/D_06_results/pad3_MAF20NA10_999thr_kmat1.csv")
write.csv(myGEMMA.npr1, "data/B05_GEMMA_les/D_06_results/npr1_MAF20NA10_999thr_kmat1.csv")
write.csv(myGEMMA.coi1, "data/B05_GEMMA_les/D_06_results/coi1_MAF20NA10_999thr_kmat1.csv")
write.csv(myGEMMA.col0, "data/B05_GEMMA_les/D_06_results/Col0_MAF20NA10_999thr_kmat1.csv")

