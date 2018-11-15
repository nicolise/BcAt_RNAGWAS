#Nicole E Soltis
#08/06/18
#---------------------------------------------------------------------
rm(list=ls())

#At-lsm
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_At/")
myPhenos <- NULL
nameddat <- NA
mydat <- NA
myPhenos <- read.table("C03_runGEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)]
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
##change these out to annotate all SNP summaries per GEMMA run
mydat <- read.table("C06_GEMMAsumm/GEMMA_top1SNPsample.txt", sep=",", row.names=NULL)
#, row.names=NULL
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
##match name here
write.csv(nameddat, "C06_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")

#--------------------------------------------------------------------------
#At-col0, coi1, npr1
##still need to do here!
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At")
myPhenos <- NULL
nameddat <- NA
mydat <- NA
##check correct phenotype
myPhenos <- read.table("col0/02_GEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)]
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
##change these out to annotate all SNP summaries per GEMMA run
<<<<<<< HEAD
mydat <- read.table("06_GEMMAsumm/npr1_GEMMA_top1SNPsample.txt", sep=",", row.names=NULL)
=======
mydat <- read.table("05_GEMMAsumm/col0_GEMMA_top100SNPsample.txt", sep=",", row.names=NULL)
>>>>>>> 97074fb931b77f5856de8808ccfd2d3ea08f4211
#, row.names=NULL
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
##match name here
<<<<<<< HEAD
write.csv(nameddat, "06_GEMMAsumm/GeneNames/npr1_GEMMA_top1SNPsample.csv")
=======
write.csv(nameddat, "05_GEMMAsumm/GeneNames/col0_GEMMA_top100SNPsample.csv")
>>>>>>> 97074fb931b77f5856de8808ccfd2d3ea08f4211
