#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data")
#for laptop setwd("~/Projects/BcSolGWAS/data/genome")
SNPs <- read.csv("02_csvPrep/hp_bin_trueMAF20_20NA.csv", row.names = 1)

Phenos <- read.csv("BOTphenotypes.csv")

#change names from genotype file to match phenotype file
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

## now only keep genotypes and phenotypes that match

#only keep phenotype rows that match SNP names
SNPMt <- as.data.frame(names(SNPs_rename))
PhenoMatch <- Phenos
PhenoMatch <- PhenoMatch[PhenoMatch$Igeno %in% SNPMt$"names(SNPs_rename)", ]

#only keep SNP rows that match phenotype names
PhenoMt <- as.data.frame(PhenoMatch[,1])
SNPMatch <- SNPs_rename
SNPs3 <- SNPs_rename[,c(1:3)]
SNPMatch <- SNPMatch[names(SNPMatch) %in% (PhenoMt$"PhenoMatch[, 1]")]
SNPMatch <- SNPMatch[ , order(names(SNPMatch))]
SNPMatch2 <- cbind(SNPs3,SNPMatch)

#REMOVE duplicate genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
SNPMatch2 <- SNPMatch2[,-9]

#sort pheno match
PhenoMatch2 <- PhenoMatch[order(PhenoMatch$Igeno),] 

#check for matching names between SNPMatch2 and PhenoMatch2
CheckNames <- PhenoMatch2[,c(1:2)]
CheckNames$SNPIgeno <- names(SNPMatch2[,c(4:96)])
CheckNames$Igeno[!(CheckNames$Igeno %in% CheckNames$SNPIgeno)] #good
CheckNames$SNPIgeno[!(CheckNames$SNPIgeno %in% CheckNames$Igeno)] #good

#now need to remove SNP columns for which all data is zero or all data is ones
SNPMatch2[which(rowSums(abs(SNPMatch2[,c(4:96)]), na.rm=T)==0),]

myvector <- rowSums(abs(SNPMatch2[,c(4:96)]), na.rm=T)
head(sort(myvector))
tail(sort(myvector))
#none: all SNPs have variation

#and, all isolates have variation
testdf <- data.frame("zero"= integer(0), "one"= integer(0))
for (i in 4:96) {
  newrow <- table(SNPMatch2[,i])
  testdf <- rbind(testdf, newrow)
}

#save them files!
write.csv(SNPMatch2, "03_bigRRinput/Domestication/hpbinSNP_bigRR_trueMAF20_50NA.csv")
write.csv(PhenoMatch2, "03_bigRRinput/Domestication/Sl_Pheno_bigRR_trueMAF20_50NA.csv")
#------------------------------------------------------------------------------
#extra things
#miniSNPs <- as.data.frame(t(miniSNPs))
#miniPhenos <- subset(Phenos, Igeno %in% SNPs_rename[0,])
