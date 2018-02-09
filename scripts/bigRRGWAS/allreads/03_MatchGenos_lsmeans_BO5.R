#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())

setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("BO5_97_iso_small/File_key_in_Bo5bamfolder_NES.csv", header=TRUE)
SNPnames <- SNPnames[,c("Isolate","names")]
names(SNPnames)[1]<- "Isolate"

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
SNPs <- read.csv("data/allreadsGWAS/BO5.10/01_prepFiles/hp_binMAF20_20NA.csv", row.names=1)
SNPs_renamed <- SNPs

Phenos <- read.csv("data/allreadsGWAS/BO5.10/01_prepFiles/lsmeans_zscale_allreads.csv")
Phenos <- Phenos[,-c(1)]

#change names from genotype file to match phenotype file
#File SNPs_renamed has columns of isolate genotypes that I want to rename
#first, remove the .variant2 from names(SNPs_renamed)
colnames(SNPs_renamed) <- sub("\\.variant2", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant3", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant", "", colnames(SNPs_renamed))
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names(SNPs_renamed) <- SNPnames[match(names(SNPs_renamed),SNPnames[,"names"]),"Isolate"] 

## now only keep genotypes and phenotypes that match
#only keep phenotype rows that match SNP names
IsosfromSNP <- as.data.frame(names(SNPs_renamed))
names(IsosfromSNP)[1] <- "IsoList"
matchedPhenos <- Phenos
intersect(matchedPhenos$Isolate, IsosfromSNP$IsoList)
setdiff(matchedPhenos$Isolate, IsosfromSNP$IsoList) #only 01.02.13
sort(unique(IsosfromSNP$IsoList))
matchedPhenos <- matchedPhenos[matchedPhenos$Isolate %in% IsosfromSNP$IsoList, ]
#yay, only drops 01.02.13

#only keep SNP rows that match phenotype names
IsosfromPheno <- as.data.frame(matchedPhenos[,1])
matchedSNPs <- SNPs_renamed
SNPs3 <- SNPs[,c(1:2)]
matchedSNPs <- matchedSNPs[names(matchedSNPs) %in% (matchedPhenos$"Isolate")] #keep 96, good
matchedSNPs <- matchedSNPs[, order(names(matchedSNPs))]
matchedSNPs2 <- cbind(SNPs3, matchedSNPs)

#REMOVE duplicate genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
matchedSNPs2 <- matchedSNPs2[,-8]

#sort pheno match
matchedPhenos2 <- matchedPhenos[order(matchedPhenos$Isolate),] 

#check for matching names between SNPMatch2 and PhenoMatch2
setdiff((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(3:97)])) #none, good
intersect((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(3:97)])) #good

#save them files!
write.csv(matchedSNPs2, "data/allreadsGWAS/BO5.10/02_bigRR/lsmeans_zscale_SNPS_MAF20.csv")
write.csv(matchedPhenos2, "data/allreadsGWAS/BO5.10/02_bigRR/lsmeans_zscale_phenos_MAF20.csv")