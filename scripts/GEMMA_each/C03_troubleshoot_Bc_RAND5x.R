#Nicole E Soltis

#confirmation: FAM file for rand2 is identical to rand1. urgh.

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS")
randFAM2 <- read.table("GEMMA_eachAt_Bc/col0/02_GEMMA_RAND02/binMAF20NA10.fam")

setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
randFAM1 <- read.table("data/GEMMA_eachAt_Bc/col0/RAND/binMAF20NA10.fam")

hist(as.numeric(as.character(randFAM2$V7)))
hist(as.numeric(as.character(randFAM1$V7)))
