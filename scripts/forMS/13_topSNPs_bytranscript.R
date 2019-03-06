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

#------------------------------------------------------------------------
#look for extreme p-values as evidence for cis
rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
top100B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
hist(-log10(top100B$p_score))
min(top100B$p_score)
hi100B <- top100B[top100B$p_score < 1e-6,]
hist(hi100B$p_score)
hist(-log10(hi100B$p_score))

top1B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt")
hist(-log10(top1B$p_score))
min(top1B$p_score)
hi1B <- top1B[top1B$p_score < 1e-6,]
hist(hi1B$p_score)
hist(-log10(hi1B$p_score))
