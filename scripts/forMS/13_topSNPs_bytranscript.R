#Nicole E Soltis
#03/01/19

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#Bc transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/04_GEMMAout/
#kmat script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A03_GEMMA_kmat.sh
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A04_runGEMMA_kmat_3genos.sh

#Need the following:
#median/ distribution of number of sig SNPs per transcript
#look for extremely low p values = evidence for cis effects

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")

#get list of phenos with many SNP > Thr
mysnplistB <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/Bc_phenos_manySNPovrThr.csv")
mysnplistA <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/At_phenos_manySNPovrThr.csv")
totnumsnpB <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/Bc_numphenosoverThr.csv")
#-----------------------------------------------------------------------------------
#calculate median SNP > thr 

#-----------------------------------------------------------------------------------
#look for extreme p-values as evidence for cis
rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
top100B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
top1B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt")

setwd("~/Projects/BcAt_RNAGWAS/data")
top100B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
top1B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt")
top100A <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top100SNPsample_genes.csv")
top1A <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes.csv")

#top 100 Bc
hist(-log10(top100B$p_score))
min(top100B$p_score)
hi100B <- top100B[top100B$p_score < 1e-6,]
hist(hi100B$p_score)
hist(-log10(hi100B$p_score))

#top 1 Bc
hist(-log10(top1B$p_score))
min(top1B$p_score)
hi1B <- top1B[top1B$p_score < 1e-6,]
hist(hi1B$p_score)
hist(-log10(hi1B$p_score))
quantile(hi1B$p_score, c(0.00, 0.01, 0.05, 0.25, 0.5))

#top 100 At
hist(-log10(top100A$p_score))
min(top100A$p_score)
hi100A <- top100A[top100A$p_score < 1e-6,]
hist(hi100A$p_score)
hist(-log10(hi100A$p_score))

#top 1 Bc
hist(-log10(top1A$p_score))
min(top1A$p_score)
hi1A <- top1A[top1B$p_score < 1e-6,]
hist(hi1A$p_score)
hist(-log10(hi1A$p_score))
quantile(hi1A$p_score, c(0.00, 0.01, 0.05, 0.25, 0.5))

library(ggplot2)
ggplot(hi1A, aes(y=-log10(p_score), x=1)) + 
  geom_violin()

ggplot(hi1B, aes(y=-log10(p_score), x=1)) + 
  geom_violin()

boxplot(-log10(hi1A$p_score), ylim=c(0, 9))
boxplot(-log10(hi1B$p_score), ylim=c(0, 9))

boxplot(-log10(hi100A$p_score), ylim=c(0, 9))
boxplot(-log10(hi100B$p_score), ylim=c(0, 9))
#use the Tukey’s method to identify the outliers ranged above and below the 1.5*IQR
boxplot.stats(-log10(hi1B$p_score))$out
boxplot.stats(-log10(hi1A$p_score))$out
