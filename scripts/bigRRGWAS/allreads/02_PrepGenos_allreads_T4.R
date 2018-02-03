#Nicole E Soltis
#09/22/17
#convert .tab data to binary .csv for all reads
#02_PrepGenos.R
#---------------------------------------------------
rm(list=ls())
#read in tab files from BcGenome
setwd("~/Projects/BcGenome/")
setwd("~/Documents/GitRepos/BcGenome")
tab20 = read.delim("data/Suzi_033016/Haploid_SNPS_97_dp6_maf20.tab")
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/")
library(tidyr)
#convert all .tab SNP files to .csv
write.table(tab20, file="01_prepFiles/snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)


#convert each file to binary
#can do for MAF5 and MAF10 too if I want
SNPsMAF20 <- read.csv("01_prepFiles/snps_maf20.csv")
mySNPs <- SNPsMAF20
#code "." as NAs
mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="."]<-NA #this is a true NA
allSNPs<- mySNPs
#replace base with 1 if match REF
#loop through it
#automatically skips NAs
mySNPs <- allSNPs
for (i in names(mySNPs[4:100])) {
  mySNPs[i][mySNPs[i]!=mySNPs$REF] <- 1
  mySNPs[i][mySNPs[i]==mySNPs$REF] <- 0
}

#remove low MAFs -- MAF20 cutoff for all
names(mySNPs)
mySNPs$Freq <- rowSums(mySNPs[,4:100] =="1")
mySNPs$Freq.0 <- rowSums(mySNPs[,4:100] =="0")
mySNPs$MAF <- (mySNPs$Freq)/ (mySNPs$Freq + mySNPs$Freq.0)
hist(mySNPs$MAF)
#check correct range
min(mySNPs$MAF, na.rm=T)
max(mySNPs$MAF, na.rm=T)
#mySNPs <- mySNPs[mySNPs$MAF <= 0.8,]

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
mySNPs <- mySNPs[mySNPs$NAcount <= 20,]
unique(mySNPs$NAcount)

#recode chrom.contig.pos
names(mySNPs)[1] <- "Chrom"
mySNPs$Chrom <- gsub(pattern = "Chromosome1$", replacement = "Chromosome1.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome2$", replacement = "Chromosome2.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome3$", replacement = "Chromosome3.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome4$", replacement = "Chromosome4.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome5$", replacement = "Chromosome5.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome6$", replacement = "Chromosome6.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome7$", replacement = "Chromosome7.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome8$", replacement = "Chromosome8.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome9$", replacement = "Chromosome9.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome10$", replacement = "Chromosome10.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome11$", replacement = "Chromosome11.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome12$", replacement = "Chromosome12.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome13$", replacement = "Chromosome13.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome14$", replacement = "Chromosome14.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome15$", replacement = "Chromosome15.0", mySNPs$Chrom)
mySNPs$Chrom <- gsub(pattern = "Chromosome16$", replacement = "Chromosome16.0", mySNPs$Chrom)
unique(mySNPs$Chrom)
mySNPs <- separate(mySNPs, Chrom, into = c("Chromosome", "Contig") )
#double check
unique(mySNPs$Chromosome)
unique(mySNPs$Contig)

#new to edit
#Reformat Chromosomes and Positions
mySNPs$Chromosome <- gsub("Chromosome", "", mySNPs$Chromosome)
mySNPs$Chromosome <- as.numeric(as.character(mySNPs$Chromosome))
mySNPs$Contig <- as.numeric(as.character(mySNPs$Contig))
names(mySNPs)[3] <- "Pos"
mySNPs$Pos <- as.numeric(as.character(mySNPs$Pos))
#remove REF and calculation variables
mySNPs <- mySNPs[,c(1:3,5:101)]

write.csv(mySNPs, "01_prepFiles/hp_binMAF20_20NA.csv")
