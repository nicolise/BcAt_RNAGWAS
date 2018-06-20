#Nicole E Soltis
#06/18/18
#AFTER  01.01.12 has been removed, are other isolates still missing phenotypes?
#and are the same SNPs dropped for each?
#----------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")

myGEMMA_1 <- read.table("data/B05_GEMMA_les/E_04_GEMMAout/binMAF20NA10_norand_kmat1_pheno1.assoc.txt")
names(myGEMMA_1) <- as.character(unlist(myGEMMA_1[1,]))
myGEMMA_1 <- myGEMMA_1[-c(1),]
myGEMMA_2 <- read.table("data/B05_GEMMA_les/E_04_GEMMAout/binMAF20NA10_norand_kmat1_pheno2.assoc.txt")
names(myGEMMA_2) <- as.character(unlist(myGEMMA_2[1,]))
myGEMMA_2 <- myGEMMA_2[-c(1),]
myGEMMA_3 <- read.table("data/B05_GEMMA_les/E_04_GEMMAout/binMAF20NA10_norand_kmat1_pheno3.assoc.txt")
names(myGEMMA_3) <- as.character(unlist(myGEMMA_3[1,]))
myGEMMA_3 <- myGEMMA_3[-c(1),]
ogSNPs <- read.csv("data/B05_GEMMA/01_PLINK/OriginalSNPdata.csv")

myGEMMA_1$chr.pos <- paste(myGEMMA_1$chr, myGEMMA_1$ps, sep=".")
myGEMMA_2$chr.pos <- paste(myGEMMA_2$chr, myGEMMA_2$ps, sep=".")
myGEMMA_3$chr.pos <- paste(myGEMMA_3$chr, myGEMMA_3$ps, sep=".")
diffSNPs <- myGEMMA_1[!(myGEMMA_1$chr.pos %in% myGEMMA_2$chr.pos),]
#all dropped SNPs identical between phenotype 1 and 2
diffSNPs <- myGEMMA_3[!(myGEMMA_3$chr.pos %in% myGEMMA_2$chr.pos),]
#and 3 vs. 1 are the same

#are there other missing phenotypes?
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
Phenos <- read.csv("LSMeanCamLes4Map_FIN.csv")
#col2 = V2 = Isolate
names(Phenos)[2] <- "V2"
#keep only lesion phenos since that's what we want to run
Phenos <- Phenos[,c("V2","Col0.Les","coi1.Les","npr1.Les","pad3.Les","tga3.Les","anac055.Les")]

sum(is.na(Phenos[,7])) #no NAs in rest of df
