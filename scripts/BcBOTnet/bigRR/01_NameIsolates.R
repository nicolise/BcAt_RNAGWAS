#092217
#Nicole E Soltis
#01_NameIsolates.R
#short script to match phenotypes of Isolates with IDs 1-97 with actual Isolate names for matching to SNP data.

#-----------------------------------------------------------------
#load data
setwd("~/Projects/BcAt_RNAGWAS/data")
IsoNames <- read.csv("IsolateKey_Vivian.csv")
BotReads <- read.csv("BotryGenes/BcBOT_result.lsm.csv")

#each plant x transcript is a unique phenotype
BotReads.w <- reshape(BotReads, idvar = "Isolate", timevar = "HostGenotype", direction = "wide")

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
BotReads.f <- merge(BotReads.w, IsoNames, by="Isolate")

#clean up columns
BotReads.f <- BotReads.f[,c(24,2:22)]
names(BotReads.f)[1] <- "Isolate"
write.csv(BotReads.f, "BOTphenotypes.csv")
