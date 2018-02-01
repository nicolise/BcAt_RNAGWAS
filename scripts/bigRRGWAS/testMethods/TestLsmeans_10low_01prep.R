#092217
#Nicole E Soltis
#testMethods/TestLsmeans_10low.R
#short script to match phenotypes of Isolates with IDs 1-97 with actual Isolate names for matching to SNP data.
#and test different scaling approaches for low lsmeans transcripts

#-----------------------------------------------------------------
rm(list=ls())
#load data
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
setwd("~/Projects/BcAt_RNAGWAS/data")
#lsm.for.GWAS comes from allreads/01_PrepPhenos_lsmeans.R
lsm.for.GWAS <- read.csv("allreadsGWAS/01_prepFiles/lsmeans_allreads.csv")

#keep only 10ish transcripts with 10% or more isolates with lsmeans < -20
#10% of 96 is 10 reads < -20
iso.names <- lsm.for.GWAS[,2]
all.my.reads <- lsm.for.GWAS[,3:9269]
my.min <- all.my.reads[1,]
for (i in (1:9267)){
  for (j in (1:10)){
    my.min[j,i] <- sort(all.my.reads[,i],partial=j)[j]   
  }
}
my.min.lowest <- my.min[,my.min[5,]<(-25)] #I'll just use this low group of 12 transcripts

#then: map these 12 transcripts 3 ways:
#a) unadjusted lsmeans
#b) lsmeans < -5 set equal to -5
#c) lsmeans z-scaled
my.names <- names(my.min.lowest)
my.low.reads <- lsm.for.GWAS[,names(lsm.for.GWAS) %in% my.names]
my.low.reads$Isolate <- lsm.for.GWAS[,2]
my.low.reads <- my.low.reads[,c(13,1:12)]

my.low.reads.cutoff <- my.low.reads
for (y in (2:13)){
  my.low.reads.cutoff[,y] <- ifelse(my.low.reads.cutoff[,y] < (-5), (-5), my.low.reads.cutoff[,y])
  names(my.low.reads.cutoff)[y] <- paste(names(my.low.reads.cutoff)[y], "cutoff", sep=".")
}

my.low.reads.z <- my.low.reads
for (z in c(2:13)){
  my.low.reads.z[,z] <- scale(my.low.reads.z[,z], center = TRUE, scale = TRUE)
  names(my.low.reads.z)[z] <- paste(names(my.low.reads.z)[z], "zscale", sep=".")
}

#now merge together
my.low.reads.all <- merge(my.low.reads, my.low.reads.cutoff, by="Isolate")
my.low.reads.all <- merge(my.low.reads.all, my.low.reads.z, by="Isolate")

write.csv(my.low.reads.all, "allreadsGWAS/01_prepFiles/lsmeans_testlowreads.csv")
my.low.reads.all <- read.csv("allreadsGWAS/01_prepFiles/lsmeans_testlowreads.csv")
#---------------------------------------------------------------------
#now prep files for bigRR
setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
names(SNPnames)[1]<- "Isolate"

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
setwd("~/Projects/BcAt_RNAGWAS/data/")
SNPs <- read.csv("allreadsGWAS/01_prepFiles/hp_binMAF20_20NA.csv", row.names = 1)
SNPs_renamed <- SNPs

Phenos <- read.csv("allreadsGWAS/01_prepFiles/lsmeans_testlowreads.csv")
Phenos <- Phenos[,-c(1)]

#change names from genotype file to match phenotype file
#File SNPs_renamed has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_renamed ) <- SNPnames[ match( names( SNPs_renamed ) , SNPnames[ , 'SNPname' ] ) ,'Isolate' ] 

## now only keep genotypes and phenotypes that match
#only keep phenotype rows that match SNP names
IsosfromSNP <- as.data.frame(names(SNPs_renamed))
names(IsosfromSNP)[1] <- "IsoList"
matchedPhenos <- Phenos
intersect(matchedPhenos$Isolate, IsosfromSNP$IsoList)
setdiff(matchedPhenos$Isolate, IsosfromSNP$IsoList)
sort(unique(IsosfromSNP$IsoList))
#have no SNPs for 1.02.13, that drops
#and rename matchedPhenos$Isolate MEAPGG to MEAP6G
levels(matchedPhenos$Isolate) <- c(levels(matchedPhenos$Isolate),"MEAP6G")
matchedPhenos$Isolate[matchedPhenos$Isolate=="MEAPGG"] <- "MEAP6G"
matchedPhenos <- matchedPhenos[matchedPhenos$Isolate %in% IsosfromSNP$IsoList, ]
#yay, only drops 01.02.13


#only keep SNP rows that match phenotype names
IsosfromPheno <- as.data.frame(matchedPhenos[,1])
matchedSNPs <- SNPs_renamed
SNPs3 <- SNPs[,c(1:3)]
matchedSNPs <- matchedSNPs[names(matchedSNPs) %in% (matchedPhenos$"Isolate")]
matchedSNPs <- matchedSNPs[, order(names(matchedSNPs))]
matchedSNPs2 <- cbind(SNPs3, matchedSNPs)

#REMOVE duplicate genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
matchedSNPs2 <- matchedSNPs2[,-9]

#sort pheno match
matchedPhenos2 <- matchedPhenos[order(matchedPhenos$Isolate),] 

#check for matching names between SNPMatch2 and PhenoMatch2
setdiff((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(4:98)]))
intersect((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(4:98)])) #good

#save them files!
write.csv(matchedSNPs2, "allreadsGWAS/02_bigRR/lsmeans_testlowreads_SNPS_MAF20.csv")
write.csv(matchedPhenos2, "allreadsGWAS/02_bigRR/lsmeans_testlowreads_phenos_MAF20.csv")


#then: evaluate outputs 2 ways
#a) correlate estimates across SNPs
#b) visual analysis of Manhattan plots across each

#-----------------------------------------------------------------------
#test for low lsmeans values
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
