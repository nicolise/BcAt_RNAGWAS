#Nicole E Soltis
#06/08/18
#D08_checkdroppedSNPs.R

#-----------------------------------------------------------------------------
#for Suzi, check MAF and missingness for SNPs that were dropped
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
myGEMMArun <- read.table("B05_GEMMA_les/D_04_ogphenos/binMAF20NA10_norand_kmat1_pheno1_Col0.assoc.Les.txt")
library("trio")
mySNPs <- read.csv("B05_GEMMA_les/D_01_PLINK/dpbinMAf20NA10ped.csv")
#full SNP list


