#Nicole E Soltis
#07/26/18
#---------------------------------------------------------------------
#B05_GEMMA_Bc done
#match GEMMA phenotypes to Gene ID
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
myPhenos <- read.table("02_GEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)] 
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
nameddat <- NA
mydat <- NA
#mydat <- read.table("05_GEMMAsumm/RawOuts/GEMMA_top100SNPsample.txt", sep=",")
mydat <- read.csv("05_GEMMAsumm/RawOuts/AllBcgenes_top10SNP_MAF20NA10_GEMMA_kmat1_Indexed.csv")
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
write.csv(nameddat, "05_GEMMAsumm/GeneNames/AllBcgenes_top10SNP_MAF20NA10_GEMMA_kmat1_Indexed.csv")

#----------------------------------------------------------------------
#B_permut
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc/B_permut/")
myPhenos <- read.table("02_GEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)]
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
nameddat <- NA
mydat <- NA

setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/B05_GEMMA_Bc/Bc_permut/")
#mydat <- read.table("05_GEMMAsumm/RawOuts/GEMMA_topSNPsample_zscale.txt", sep=",", row.names=NULL)
##change these out to annotate all SNP summaries per GEMMA run
mydat <- read.table("05_GEMMAsumm/RawOuts/GEMMA_topSNPsample_zscale.txt", sep=",")
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
write.csv(nameddat, "05_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")
#--------------------------------------------------------------------
#Bc_col0 / coi1 / npr1 
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc")
myPhenos <- NULL
nameddat <- NA
mydat <- NA
myPhenos <- read.table("npr1/02_GEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)]
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
##change these out to annotate all SNP summaries per GEMMA run
mydat <- read.table("06_GEMMAsumm/npr1_GEMMA_top1SNPsample.txt", sep=",", row.names=NULL)
#, row.names=NULL
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
##match name here
write.csv(nameddat, "06_GEMMAsumm/GeneNames/npr1_GEMMA_top1SNPsample.csv")
#---------------------------------------------------------------
