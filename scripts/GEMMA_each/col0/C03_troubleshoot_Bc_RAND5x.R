#Nicole E Soltis

#confirmation: FAM file for rand2 is identical to rand1. urgh.

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS")
randFAM2 <- read.table("GEMMA_eachAt_Bc/col0/02_GEMMA_RAND02/binMAF20NA10.fam")

setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
randFAM1 <- read.table("data/GEMMA_eachAt_Bc/col0/RAND/binMAF20NA10.fam")

#based on hist they're identical
hist(as.numeric(as.character(randFAM2$V7)))
hist(as.numeric(as.character(randFAM1$V7)))

#but different order?
plot(as.numeric(as.character(randFAM2$V8)), as.numeric(as.character(randFAM1$V8)))
#not identical here

#what about the outputs? -- yes, identical. RERUNNING rand2 with unique phenotype order.
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS")
rand2assc1 <- read.table("GEMMA_eachAt_Bc/04_GEMMAout/col0_5Rand/rand2/col0_MAF20NA10_RAND02_1.assoc.txt", header = TRUE)
rand1assc1 <- read.table("GEMMA_eachAt_Bc/04_GEMMAout/col0_5Rand/rand1/col0_MAF20NA10_RAND_1.assoc.txt", header = TRUE)
names(rand1assc1)
hist(rand1assc1$p_score)
hist(rand2assc1$p_score)
