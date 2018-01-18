#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())

setwd("~/Projects/BcGenome/data")
SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
names(SNPnames)[1]<- "Isolate"

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
setwd("~/Projects/BcAt_RNAGWAS/data/")
SNPs <- read.csv("allreadsGWAS/01_prepFiles/hp_binMAF20_20NA.csv", row.names = 1)
#SNPs <- read.csv("02_MatchGenos/hp_binMAF10_20NA.csv", row.names = 1)
#SNPs <- read.csv("02_MatchGenos/hp_binMAF05_20NA.csv", row.names = 1)
SNPs_renamed <- SNPs

Phenos <- read.csv("allreadsGWAS/01_prepFiles/col0_allreads.csv")
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
write.csv(matchedSNPs2, "allreadsGWAS/02_bigRR/allreads_SNPS_MAF20.csv")
write.csv(matchedPhenos2, "allreadsGWAS/02_bigRR/allreads_phenos_MAF20.csv")
