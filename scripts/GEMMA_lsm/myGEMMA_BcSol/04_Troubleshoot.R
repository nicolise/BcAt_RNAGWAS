#04_troubleshoot

#troubleshoot missingness
#--------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

#this fails, not sure why
#myPed <- read.delim("GEMMA_files/02_csvPrep/fulldata/01_PLINK_55chrom/dpcharMAF20NA10.ped", sep=" ")
myPed <- read.csv("GEMMA_files/02_csvPrep/fulldata/PED_dpcharMAF20NA10.csv")
names(myPed)[1:7]
head(myPed[,1:30], 12)
#NA is coded as <NA>
is.na(myPed[7,26])
#Individuals are rows, SNPs are columns (but each 2 columns = 1 SNP)

myMAP <- 

#This will remove all SNP columns containing at least one NA:
myPed2 <- myPed[ , colSums(is.na(myPed)) == 0]
#myPed2 has only SNPs with complete data for all isolates
#just double-check that all NAs are gone:
myPed2.2 <- myPed2[rowSums(is.na(myPed2))==0,]

myPed3 <- myPed[ rowSums(is.na(myPed)) == 0, ]
#there are no isolates with data for all SNPs

#why not try to add in phenotypes at this step?
#I can only add one phenotype column before going ped -> bed
#so for now I will only add the first phenotype (LA1547) as a trial run
myPhenos <- read.csv("GEMMA_files/02_csvPrep/Sl_Pheno_bigRR_trueMAF20_10NA.csv")
myKey <- read.csv("GEMMA_files/02_csvPrep/Key_SNPnames.csv")
myKey <- myKey[-c(1:3),]
myKey <- myKey[,c(2,4)]
names(myKey)[1] <- paste("Igeno")
names(myFam)[2] <- paste("SNPname")

#first: add column with isolate nicknames from Suzi to phenos
myPhenos.rn <- merge(myPhenos, myKey, by="Igeno")
myPhenos.rn <- myPhenos.rn[,c(15,3)]
names(myPhenos.rn)[1] <- "Isolate"

#then: merge myPhenos file to myPed
#for now only LA1547 phenotype is in use
names(myPed)[1:7]
head(myPed2.2.p[,1:30], 12)
myPed2.2.p <- merge(myPed2.2, myPhenos.rn, by="Isolate")
#reorder to make this sensible
myPed.fin <- myPed2.2.p[,c(2:6,504112,8:504111)]
head(myPed.fin[,1:30], 12)
#add all 12 phenotypes. Each time I run GEMMA, specify which column (col 6 = 1, col 7 = 2, etc.)

library(caroline)
write.delim(myPed.fin, "GEMMA_files/02_csvPrep/noMissing/dpcharMAF20NA10.ped", quote = FALSE, col.names = F, row.names = FALSE, sep = "\t")


write.csv(myPed2.2, "GEMMA_files/02_csvPrep/noMissing/PED_dpcharMAF20NA10.csv")
