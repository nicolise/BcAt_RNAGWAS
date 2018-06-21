#092217
#Nicole E Soltis
#01_NameIsolates.R
#short script to match phenotypes of Isolates with IDs 1-97 with actual Isolate names for matching to SNP data.

#-----------------------------------------------------------------
rm(list=ls())
library(lsmeans); library(data.table)
#load data
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
IsoNames <- read.csv("Vivian_Bc/IsolateKey_Vivian.csv")
MyReads <- read.csv("Vivian_At/data from Wei 20160321/S2_ModelLsMeanAdjusted.csv")
MyReads <- MyReads[,-c(1)]

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(23960,2,3:23959)]
names(MyReads)[1] <- "Isolate"

#approach: new lsmeans including the effect of HostGenotype
full.lsm.outputs=NULL

#lsmeans::lsmeans is being deprecated, now using emmeans
library(emmeans)

#first run once to get the right isolate names because I'm lazy
#fix this later
for (i in c(4)) {
  #Infected/Isolate
  MyReads.lm <- lm(MyReads[,i] ~ Infected/Isolate + HostGenotype, data=MyReads)
  #MyReads.lsm <- as.data.frame(print(lsmeans::lsmeans(MyReads.lm, "Isolate")))
  MyReads.lsm <- as.data.frame(print(emmeans::emmeans(MyReads.lm, "Isolate")))
  df <- as.data.frame(print(MyReads.lsm))
  #setDT(df, keep.rownames = T)[]
  df$transcript <- names(MyReads)[i]
  full.lsm.outputs = rbind(full.lsm.outputs, df)
  lsm.for.GWAS = (df)
  names(lsm.for.GWAS)[i-1] <- names(MyReads)[i]
}
lsm.for.GWAS= as.data.frame(MyReads.lsm[,1])
names(lsm.for.GWAS)[1] <- "Isolate"
#run for 4:23959
mytime1 <- (Sys.time())
#takes 5 mins for 1000 -> took 11 hours for 23957
for (i in c(4:length(MyReads))) {
  MyReads.lm <- lm(MyReads[,i] ~ Infected/Isolate + HostGenotype, data=MyReads)
  #MyReads.lsm <- as.data.frame(print(lsmeans::lsmeans(MyReads.lm, "Isolate")))
  MyReads.lsm <- as.data.frame(print(emmeans::emmeans(MyReads.lm, "Isolate")))
  df <- as.data.frame(print(MyReads.lsm))
  #setDT(df, keep.rownames = T)[]
  df$transcript <- names(MyReads)[i]
  full.lsm.outputs = rbind(full.lsm.outputs, df)
  lsm.for.GWAS = cbind(lsm.for.GWAS,df$emmean)
  names(lsm.for.GWAS)[i-2] <- names(MyReads)[i]
}
print(mytime1)
print(Sys.time())

#try z-scaling 
mytime1 <- Sys.time()
my.lsm.z <- lsm.for.GWAS[,-c(1)]
for (i in c(1:(length(my.lsm.z)))){
  my.lsm.z[,i] <- scale(my.lsm.z[,i], center = TRUE, scale = TRUE)
}
my.lsm.z$Isolate <- lsm.for.GWAS[,1]
print(mytime1)
print(Sys.time())
my.lsm.z <- my.lsm.z[,c(9268,1:9267)]
  
#write out file
write.csv(full.lsm.outputs, "B05_GEMMA/C01_AtPrepFiles/ANOVA_lsmoutputs_allreads_At.csv")
write.csv(lsm.for.GWAS, "B05_GEMMA/C01_AtPrepFiles/lsmeans_allreads_At.csv")
write.csv(my.lsm.z, "B05_GEMMA/C01_AtPrepFiles/lsmeans_zscale_allreads_At.csv")

#-----------------------------------------------------------------------
## to customize for this script
#test for low lsmeans values
my.lsm <- read.csv("B05_GEMMA/C01_AtPrepFiles/lsmeans_allreads_At.csv")
names(my.lsm)[1:10]
min(my.lsm)
all.my.reads <- my.lsm[,3:length(my.lsm)]
my.mins <- apply(all.my.reads, 2, min)
hist(my.mins)
min(my.mins) #down to -28
my.maxs <- apply(all.my.reads, 2, max)
hist(my.maxs)
max(my.maxs) #up to 12
sum( my.mins < 0)
sum( my.mins < -1)
sum( my.mins < -2)
sum( my.mins < -5)
my.mins <- as.data.frame(my.mins)
my.maxs <- as.data.frame(my.maxs)
library(ggplot2)
ggplot(data=my.mins, aes(my.mins)) + geom_histogram(bins=40)
ggplot(data=my.maxs, aes(my.maxs)) + geom_histogram(bins=40)

my.2.min <- all.my.reads[1,]
for (i in (1:9267)){
  my.2.min[,i] <- sort(all.my.reads[,i],partial=2)[2]
}
my.2nd.min <- (as.numeric(my.2.min[1,]))
my.2nd.min <- as.data.frame(my.2nd.min)
ggplot(data=my.2nd.min, aes(my.2nd.min)) + geom_histogram(bins=40)

my.6.min <- all.my.reads[1,]
for (i in (1:9267)){
  my.6.min[,i] <- sort(all.my.reads[,i],partial=6)[6]
}
my.6th.min <- (as.numeric(my.6.min[1,]))
my.6th.min <- as.data.frame(my.6th.min)
ggplot(data=my.6th.min, aes(my.6th.min)) + geom_histogram(bins=40)

#try z-scaling 
my.lsm.z <- my.lsm[,-c(1)]
for (i in c(2:9268)){
  my.lsm.z[,i] <- scale(my.lsm.z[,i], center = TRUE, scale = TRUE)
  #names(my.lsm.z)[i] <- print(names(my.lsm)[i])
}
my.mins.z <- apply(my.lsm.z, 2, min)
my.mins.z <- as.numeric(my.mins.z[2:9268])
hist(my.mins.z)
min(my.mins.z) 
