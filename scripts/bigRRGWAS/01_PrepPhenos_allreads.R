#092217
#Nicole E Soltis
#01_NameIsolates.R
#short script to match phenotypes of Isolates with IDs 1-97 with actual Isolate names for matching to SNP data.

#-----------------------------------------------------------------
rm(list=ls())
#load data
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
setwd("~/Projects/BcAt_RNAGWAS/data")
IsoNames <- read.csv("IsolateKey_Vivian.csv")
MyReads <- read.csv("Vivian_Bc/result.lsm.csv")
MyReads <- MyReads[,-c(1)]

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(9270,2,3:9269)]
names(MyReads)[1] <- "Isolate"

#each plant x transcript is a unique phenotype
#split by host
splitlist <- split(MyReads, MyReads$HostGenotype)
unique(MyReads$HostGenotype)
MyReads.coi1 <- splitlist[[1]]
MyReads.col0 <- splitlist[[2]]
MyReads.npr1 <- splitlist[[3]]

#other approach: new lsmeans including the effect of HostGenotype

full.lsm.outputs=NULL
lsm.for.GWAS= as.data.frame(MyReads.lsm[,1])
names(lsm.for.GWAS)[1] <- "Isolate"
library(data.table)
#run for 3:9269
mytime1 <- (Sys.time())
#takes 5 mins for 1000 -> predict 45 mins for 9k
for (i in c(3:9269)) {
  MyReads.lm <- lm(MyReads[,i] ~ Isolate + HostGenotype, data=MyReads)
  MyReads.lsm <- as.data.frame(print(lsmeans(MyReads.lm, "Isolate")))
  df <- as.data.frame(print(MyReads.lsm))
  #setDT(df, keep.rownames = T)[]
  df$transcript <- names(MyReads)[i]
  full.lsm.outputs = rbind(full.lsm.outputs, df)
  lsm.for.GWAS = cbind(lsm.for.GWAS,df$lsmean)
  names(lsm.for.GWAS)[i-1] <- names(MyReads)[i]
}
print(mytime1)
print(Sys.time())

#write out file
write.csv(MyReads.col0, "allreadsGWAS/01_prepFiles/col0_allreads.csv")
write.csv(full.lsm.outputs, "allreadsGWAS/01_prepFiles/ANOVA_lsmoutputs_allreads.csv")
write.csv(lsm.for.GWAS, "allreadsGWAS/01_prepFiles/lsmeans_allreads.csv")
