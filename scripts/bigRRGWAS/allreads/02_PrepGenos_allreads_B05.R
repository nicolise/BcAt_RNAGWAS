#Nicole E Soltis
#09/22/17
#convert .tab data to binary .csv for all reads
#02_PrepGenos.R
#---------------------------------------------------
rm(list=ls())
#read in tab files from BcGenome
setwd("~/Documents/GitRepos/BcGenome")
#minor allele count (MAC) > 20%
tab20 = read.delim("data/BO5_97_iso_small/snps_mac20_BO5_97.tab")
#could also read in haploidized vcf
#data/BO5_97_iso_small/mac20_hap_97_BO5.vcf
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
library(tidyr)
#convert all .tab SNP files to .csv
write.table(tab20, file="data/allreads_bigRR/BO5.10/01_prepFiles/snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)

#convert each file to binary
SNPsMAF20 <- read.csv("data/allreads_bigRR/B05.10/01_prepFiles/snps_maf20.csv")
mySNPs <- SNPsMAF20
#code "." as NAs
# str(mySNPs)
#make these characters instead of factors
mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="./"] <- NA #this is a true NA
mySNPs[mySNPs=="G/"] <- "G"
mySNPs[mySNPs=="C/"] <- "C"
mySNPs[mySNPs=="A/"] <- "A"
mySNPs[mySNPs=="T/"] <- "T"
allSNPs<- mySNPs

#replace base with 1 if match REF
#loop through it
#automatically skips NAs
mySNPs <- allSNPs
for (i in names(mySNPs[4:100])) {
  mySNPs[i][mySNPs[i]!=mySNPs$REF] <- 1
  mySNPs[i][mySNPs[i]==mySNPs$REF] <- 0
}

#check for low MAFs -- MAF20 cutoff for all
#kept all
names(mySNPs)
mySNPs$Freq <- rowSums(mySNPs[,4:100] =="1")
mySNPs$Freq.0 <- rowSums(mySNPs[,4:100] =="0")
mySNPs$MAF <- (mySNPs$Freq)/ (mySNPs$Freq + mySNPs$Freq.0)
hist(mySNPs$MAF)
#check correct range
min(mySNPs$MAF, na.rm=T)
max(mySNPs$MAF, na.rm=T)
#mySNPs <- mySNPs[mySNPs$MAF <= 0.8,]

#and check for missing SNP calls
mySNPs$Freq.1 <- rowSums(mySNPs[,4:100] =="1", na.rm=T)
mySNPs$Freq.0 <- rowSums(mySNPs[,4:100] =="0", na.rm=T)
mySNPs$NAcount <- 97 - (mySNPs$Freq.1 + mySNPs$Freq.0)
mySNPs$Freq <- (mySNPs$Freq.1)/ (mySNPs$Freq.1 + mySNPs$Freq.0)
hist(mySNPs$Freq)
hist(mySNPs$NAcount)

#now, make choices:
#omit loci with low info?
#total of 97 isolates. A SNP with data in 90% of isolates is present in:
97*.9 #>88 isolates
#meaning NA in:
97*.1 #<10 isolates
#or for data in 80% of isolates:
97*.2 #<20 isolates
#data in 50% of isolates:
97*.5 #<49 isolates

#ok, remove all SNPs with NA in > 20 isolates
mySNPs <- mySNPs[mySNPs$NAcount <= 20,]
unique(mySNPs$NAcount)

#no contigs for BO5.10, cool!
names(mySNPs)[1] <- "Chrom"

#new to edit
#Reformat Chromosomes and Positions
mySNPs$Chrom <- gsub("chr", "", mySNPs$Chrom)
mySNPs$Chrom <- as.numeric(as.character(mySNPs$Chrom))
names(mySNPs)[2] <- "Pos"
mySNPs$Pos <- as.numeric(as.character(mySNPs$Pos))
#remove REF and calculation variables
mySNPs <- mySNPs[,c(1:2,4:100)]

#check for duplicated SNPs
checkSNPs <- mySNPs
checkSNPs$chr.pos <- paste(checkSNPs$Chrom, checkSNPs$Pos, sep=".")
checkSNPs$dup <- duplicated(checkSNPs$chr.pos)
#none duplicated

write.csv(mySNPs, "data/allreadsGWAS/BO5.10/01_prepFiles/hp_binMAF20_20NA.csv")
