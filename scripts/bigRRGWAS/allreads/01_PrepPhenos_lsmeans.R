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

#-----------------------------------------------------------------------
#test for low lsmeans values
my.lsm <- read.csv("allreadsGWAS/01_prepFiles/lsmeans_allreads.csv")
names(my.lsm)[1:10]
min(my.lsm)
all.my.reads <- my.lsm[,3:9269]
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
