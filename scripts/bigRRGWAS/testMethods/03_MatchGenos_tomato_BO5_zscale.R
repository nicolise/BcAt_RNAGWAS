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

setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files")
IPhenos <- read.csv("02_csvPrep/phenos/NewModel0711/BcSl_lsmeans_forbigRR.csv")
DPhenos <- read.csv("02_csvPrep/phenos/Domestication/BcSl_lsmeans_domest_forbigRR.csv")
Phenos <- merge(DPhenos,IPhenos, by="Igeno")

for (i in c(2:16)){
  Phenos[,i] <- scale(Phenos[,i], center = TRUE, scale = TRUE)
}
#change names from genotype file to match phenotype file
#File SNPs_renamed has columns of isolate genotypes that I want to rename
#first, remove the .variant2 from names(SNPs_renamed)
colnames(SNPs_renamed) <- sub("\\.variant2", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant3", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant", "", colnames(SNPs_renamed))
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names(SNPs_renamed) <- SNPnames[match(names(SNPs_renamed),SNPnames[,"names"]),"Isolate"] 
#rename MEAPGG in SNPs_renamed to MEAP6G
names(SNPs_renamed)[93] <- "MEAP6G"
## now only keep genotypes and phenotypes that match
#only keep phenotype rows that match SNP names
IsosfromSNP <- as.data.frame(names(SNPs_renamed))
names(IsosfromSNP)[1] <- "IsoList"
matchedPhenos <- Phenos
intersect(matchedPhenos$Igeno, IsosfromSNP$IsoList)
setdiff(matchedPhenos$Igeno, IsosfromSNP$IsoList) #only 01.02.05, 94.4 (no sequencing)

sort(unique(IsosfromSNP$IsoList))
matchedPhenos <- matchedPhenos[matchedPhenos$Igeno %in% IsosfromSNP$IsoList, ]
#yay, only drops 01.02.05, 94.4

#only keep SNP rows that match phenotype names
IsosfromPheno <- as.data.frame(matchedPhenos[,1])
matchedSNPs <- SNPs_renamed
SNPs3 <- SNPs[,c(1:2)]
matchedSNPs <- matchedSNPs[names(matchedSNPs) %in% (matchedPhenos$"Igeno")] #keep 94, good
matchedSNPs <- matchedSNPs[, order(names(matchedSNPs))]
matchedSNPs2 <- cbind(SNPs3, matchedSNPs)

#REMOVE duplicate genotype column with 01.01.06.1
#remove SNP column "X1.01.06.1"
matchedSNPs2 <- matchedSNPs2[,-8]

#sort pheno match
matchedPhenos2 <- matchedPhenos[order(matchedPhenos$Igeno),] 

#check for matching names between SNPMatch2 and PhenoMatch2
setdiff((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(3:95)])) #none, good
intersect((matchedPhenos2[,c(1)]), names(matchedSNPs2[,c(3:95)])) #all 93, good

#save them files!
setwd("~/Projects/BcAt_RNAGWAS/")
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
write.csv(matchedSNPs2, "data/testMethods/BO5tomato/02_bigRR/tomato_zscale_SNPS_MAF20.csv")
write.csv(matchedPhenos2, "data/testMethods/BO5tomato/02_bigRR/tomato_zscale_phenos_MAF20.csv")
