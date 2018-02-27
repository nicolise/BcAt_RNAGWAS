#Nicole E Soltis
#Prep .inp file for fastPHASE
#in order to call haplotypes in cis with genes
#----------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#use same SNP set from B05.10 bigRR
mySNPs <- read.csv("allreadsGWAS/BO5.10/01_prepFiles/hp_binMAF20_20NA.csv")
SNPs_renamed <- mySNPs
#change names from genotype file to match phenotype file
colnames(SNPs_renamed) <- sub("\\.variant2", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant3", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant", "", colnames(SNPs_renamed))
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("BO5_97_iso_small/File_key_in_Bo5bamfolder_NES.csv", header=TRUE)
SNPnames <- SNPnames[,c("Isolate","names")]
names(SNPnames)[1]<- "Isolate"
names(SNPs_renamed) <- SNPnames[match(names(SNPs_renamed),SNPnames[,"names"]),"Isolate"] 
mySNPs <- SNPs_renamed[,-c(1)]
names(mySNPs)[1] <- "Chrom"
names(mySNPs)[2] <- "Pos"
mySNPs <- mySNPs[,-c(38)] # remove 1.01.06.1



#line 3n: name of individual n
#line 3n + 1: list of genotypes (SNPs, in order) for individual n
  #use ? for missing info
#line 3n + 2: duplicate line above to make fake homozygous diploids

#replace all NA with ?
mySNPs[is.na(mySNPs)] <- "?"

#only keep SNPs in region of interest

#transpose data frame for region of interest
myTSNPs <- as.data.frame(t(mySNPs))

#start writing out file
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#first line: no.individuals = 96 
myline1 <- 96
write(myline1,file="fastPHASE/mypractice.txt",append=FALSE)
#second line: no.SNPsites = 271749
myline2 <- 271749
write(myline2,file="fastPHASE/mypractice.txt",append=TRUE)
myline3 <- rownames(myTSNPs)[6]
write(myline3,file="fastPHASE/mypractice.txt",append=TRUE)
myline4 <- myTSNPs[6,]
write(myline4, file="fastPHASE/mypractice.txt",append=TRUE)
# myline4 <- as.data.frame(mySNPs[,6])
# myline4 <- cbind(myline4, myline4)
# myline5 <- as.data.frame(t(myline4))
# myline6 <- as.vector(myline5[1,])
# write(myline6,file="fastPHASE/mypractice.txt",append=TRUE)
