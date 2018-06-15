#Nicole E Soltis
#06/08/18
#D08_checkdroppedSNPs.R

#-----------------------------------------------------------------------------
#for Suzi, check MAF and missingness for SNPs that were dropped
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
myGEMMArun <- read.table("B05_GEMMA_les/D_04_ogphenos/binMAF20NA10_norand_kmat1_pheno1_Col0.assoc.Les.txt")
library("trio")
#full SNP list
myMAP <- read.table("B05_GEMMA_les/D_01_PLINK/dpbinMAF20NA10.map")
ogSNPs <- read.csv("B05_GEMMA/01_PLINK/OriginalSNPdata.csv")

head(myMAP)
myMAP <- myMAP[,-c(3)]
names(myMAP) <- c("chr", "snpname", "position")

names(myGEMMArun) <- c("chr", "snpname", "position", "n_mis", "n_obs", "allele1", "allele0", "af", "beta", "se", "p_wald", "p_lrt", "p_score")
myGEMMArun <- myGEMMArun[-c(1),]

ogSNPs$chr.pos <- paste(ogSNPs$Chrom, ogSNPs$Pos, sep=".")
myGEMMArun$chr.pos <- paste(myGEMMArun$chr, myGEMMArun$position, sep=".")
dropSNPs <- ogSNPs[!(ogSNPs$chr.pos %in% myGEMMArun$chr.pos),]
#35595 SNPs dropped

dropSNPs <- dropSNPs[,-c(1)]
dropSNPs$missingness <- rowSums(is.na(dropSNPs))

dropSNPs$Freq.1 <- rowSums(dropSNPs[,3:99] =="1", na.rm=T)
dropSNPs$Freq.0 <- rowSums(dropSNPs[,3:99] =="0", na.rm=T)
#+1 for mystery individual who was dropped from GEMMA
dropSNPs$Freq <- (dropSNPs$Freq.1)/ (dropSNPs$Freq.1 + dropSNPs$Freq.0 + dropSNPs$missingness + 1)
hist(dropSNPs$Freq)

#max 9 SNPs can be missing for 96 total individuals
  #unsure which individual was dropped
dropSNPs$failed <- ifelse(dropSNPs$missingness > 9, "missingness",
                          ifelse(dropSNPs$Freq < 0.2, "lowMAF","unknown"))
                          
write.csv(dropSNPs, "B05_GEMMA_les/D_01_PLINK/DroppedSNPs_GEMMA.csv")
dropSNPs <- read.csv("data/B05_GEMMA_les/D_01_PLINK/DroppedSNPs_GEMMA.csv")

#visualize LD
#on laptop
