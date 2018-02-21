#03_FAM_addPhenotypes
#Nicole E Soltis
#072717
#add actual phenotype column to .fam file
#--------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

myFam <- read.delim("GEMMA_files/02_csvPrep/fulldata/01_PLINK_55chrom/dpcharMAF20NA10_out.fam", sep=" ")
myPhenos <- read.csv("GEMMA_files/02_csvPrep/Sl_Pheno_bigRR_trueMAF20_10NA.csv")
myKey <- read.csv("GEMMA_files/02_csvPrep/Key_SNPnames.csv")
myKey <- myKey[-c(1:3),]
myKey <- myKey[,c(2,4)]
names(myKey)[1] <- paste("Igeno")
names(myFam)[2] <- paste("SNPname")

#first: add column with isolate nicknames from Suzi to phenos
myPhenos.rn <- merge(myPhenos, myKey, by="Igeno")

#then: merge myPhenos file to myFam
myFam.p <- merge(myFam, myPhenos.rn, by="SNPname")

#add all 12 phenotypes. Each time I run GEMMA, specify which column (col 6 = 1, col 7 = 2, etc.)

myFam.p <- myFam.p[,-c(6:8)]

write.delim(myFam.p, "GEMMA_files/02_csvPrep/fulldata/PLINK_55chrom/FAM_dpcharMAF20NA10_55chrom_Phenos.fam", quote = FALSE, col.names = F, row.names = FALSE, sep = "\t")


#now, need to edit kinship file to match number of phenotypes:
