#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
SNPs <- read.csv("02_MatchGenos/hp_binMAF20_20NA.csv", row.names = 1)
#SNPs <- read.csv("02_MatchGenos/hp_binMAF10_20NA.csv", row.names = 1)
#SNPs <- read.csv("02_MatchGenos/hp_binMAF05_20NA.csv", row.names = 1)
SNPs_rename <- SNPs

SNPnames <- read.csv("02_MatchGenos/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,4)]
names(SNPnames)[1]<- "Isolate"

Phenos <- read.csv("02_MatchGenos/BOTphenotypes.csv")
Phenos <- Phenos[,-c(1)]

#change names from genotype file to match phenotype file
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'Isolate' ] 

## now only keep genotypes and phenotypes that match

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch <- Phenos
PhenoMatch <- PhenoMatch[PhenoMatch$Isolate %in% SNPMt$"names(SNPs_rename)", ]
#only lost 01.02.13 -- fine

#only keep SNP rows that match phenotype names
PhenoMt <- as.data.frame(PhenoMatch[,1])
SNPMatch <- SNPs_rename
SNPs3 <- SNPs[,c(1:3)]
SNPMatch <- SNPMatch[names(SNPMatch) %in% (PhenoMt$"PhenoMatch[, 1]")]
SNPMatch <- SNPMatch[ , order(names(SNPMatch))]
SNPMatch2 <- cbind(SNPs3,SNPMatch)

#REMOVE duplicate genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
SNPMatch2 <- SNPMatch2[,-9]

#sort pheno match
PhenoMatch2 <- PhenoMatch[order(PhenoMatch$Isolate),] 

#check for matching names between SNPMatch2 and PhenoMatch2
CheckNames <- as.data.frame(PhenoMatch2[,c(1)])
CheckNames$SNPIgeno <- names(SNPMatch2[,c(4:98)])
CheckNames$`PhenoMatch2[, c(1)]`[!(CheckNames$`PhenoMatch2[, c(1)]` %in% CheckNames$SNPIgeno)] #good
CheckNames$SNPIgeno[!(CheckNames$SNPIgeno %in% CheckNames$`PhenoMatch2[, c(1)]`)] #good


#save them files!
write.csv(SNPMatch2, "03_bigRR/BOT_genos_MAF05.csv")
write.csv(PhenoMatch2, "03_bigRR/BOT_phenos_MAF05.csv")
#------------------------------------------------------------------------------
#extra things
#miniSNPs <- as.data.frame(t(miniSNPs))
#miniPhenos <- subset(Phenos, Igeno %in% SNPs_rename[0,])
