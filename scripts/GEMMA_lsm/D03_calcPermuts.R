#Nicole E Soltis
#07/03/18

#-------------------------------------------------------------------
rm(list=ls())

#get SNP number from allmySNP = 237096
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
allmySNP <- read.table("04_GEMMAoutput/binMAF20NA10_PLINK_1.assoc.txt")

#get gene number from allmyGenes = 9270
setwd("~/Projects/BcAt_RNAGWAS/data/Vivian_Bc")
allmyGenes <- read.csv("result.lsm.csv")

alldat <- as.data.frame(1:237096)
names(alldat)[1] <- "SNP"
for (i in 1:1000){
  mylist <- sample(1:237096,9270,replace=T)
  mydat <- as.data.frame(table(mylist))
  names(mydat)[1] <- "SNP"
  names(mydat)[2] <- paste("Freq",i,sep="_")
  alldat <- merge(alldat, mydat, by = "SNP", all=TRUE)
}
max(alldat[,2:length(alldat)], na.rm=T)
max(alldat[,2:1044], na.rm=T)

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc/05_GEMMAsumm")
write.csv(alldat, "RandSNPoverlap_simulated.csv")

cbind()
names(mydat)[1] <- "SNP"
max(mydat$Freq)
