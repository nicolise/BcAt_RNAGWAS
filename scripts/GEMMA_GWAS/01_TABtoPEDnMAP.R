#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
library(tidyr)
#convert each file to binary
SNPsMAF20 <- read.csv("data/allreadsGWAS/BO5.10/01_prepFiles/snps_maf20.csv")
mySNPs <- SNPsMAF20
#make these characters instead of factors
mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="./"] <- NA #this is a true NA
mySNPs[mySNPs=="G/"] <- "G"
mySNPs[mySNPs=="C/"] <- "C"
mySNPs[mySNPs=="A/"] <- "A"
mySNPs[mySNPs=="T/"] <- "T"
allSNPs<- mySNPs

#remove all SNPs with NA in > 20 isolates to match 02_PrepGenos_allreads_B05.R
mySNPs$NAcount <- apply(mySNPs[,c(4:100)], 1, function(x) sum(is.na(x)))
mySNPs <- mySNPs[mySNPs$NAcount <= 20,]
unique(mySNPs$NAcount)

#may as well rename isolates here
setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("BO5_97_iso_small/File_key_in_Bo5bamfolder_NES.csv", header=TRUE)
SNPnames <- SNPnames[,c("Isolate","names")]
names(SNPnames)[1]<- "Isolate"
SNPs_renamed <- mySNPs
colnames(SNPs_renamed) <- sub("\\.variant2", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant3", "", colnames(SNPs_renamed))
colnames(SNPs_renamed) <- sub("\\.variant", "", colnames(SNPs_renamed))
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names(SNPs_renamed) <- SNPnames[match(names(SNPs_renamed),SNPnames[,"names"]),"Isolate"] 
names(SNPs_renamed)[1] <- "Chrom"
names(SNPs_renamed)[2] <- "Pos"
names(SNPs_renamed)[3] <- "REF"
#and now for making PED format for PLINK!
  #do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#add a first column: FAM1 (no info on isolate families)
#second column: isolate ID
#third column: father ID (a column of zeros)
#fourth column: mother ID (a column of zeros)
#fifth column: individual sex = 1 (all assumed same)
#sixth  column: binary  phenotype (all = 1)
#fix column order
mySNPs2 <- SNPs_renamed[,-c(1:3,101)]
#remove duplicate for 1.01.06.1
mySNPs2 <- mySNPs2[,-c(36)]

#for PED, NA must be replaced with 0 for genotypes, else NA will be read as an allele
mySNPs2[is.na(mySNPs2)] <- 0

#turn all SNPs to "diploid"
#haha, it takes 4 days to do this as a "for" loop (for each row, rbind twice)
#because is.na <-0 before this step, there should be NO heterozygous SNP calls
#this is super fast:
mySNPs3 <- mySNPs2[rep(1:nrow(mySNPs2),each=2),] 

#transpose and format for PED
mySNPs3 <- mySNPs3[,-c(1)]
mySNPs4 <- as.data.frame(t(mySNPs3))
#add binary phenotype = 1 (6)
mySNPs4 <- cbind("Pheno" = 1, mySNPs4)
#add individual sex = 1 (5)
mySNPs4 <- cbind("sex" = 1, mySNPs4)
#add Mother = 0 (4)
mySNPs4 <- cbind("Mother" = 0, mySNPs4)
#add Father = 0 (3)
mySNPs4 <- cbind("Father" = 0, mySNPs4)
#turn row names into column 2
mySNPs4 <- cbind(rownames(mySNPs4), mySNPs4)
colnames(mySNPs4)[1] <- 'Isolate'
#add the fam column (1)
mySNPs4 <- cbind("FAM" = "FAM1", mySNPs4)
myPED <- mySNPs4

#add a phenotype for PED? 
#NA is fine for missing phenotypes
#since many phenotypes, just add as consecutive columns to *.fam, and run GEMMA in a loop over phenotypes

#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- SNPs_renamed[,c("Chrom","Pos")]
#remove "Chromosome" from X.CHROM
myMAP2 <- as.data.frame(lapply(myMAP, function(x) {
                 gsub("chr", "", x)
              }))
myMAP2$SNPID <- paste("SNP",myMAP2$Pos, sep="")
myMAP2$SNPcM <- 0
myMAP2 <- myMAP2[,c(1,3,4,2)]
write.table(myMAP2, "data/B05_GEMMA/01_PLINK/dpcharMAF20NA10.map", row.names=FALSE, col.names=FALSE)
#add a column of "SNP identifiers" in excel and remove headers

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
write.csv(mySNPs3, "data/B05_GEMMA/01_PLINK/dp_charMAF20_10NA.csv")
write.csv(mySNPs, "data/B05_GEMMA/01_PLINK/hp_charMAF20_10NA.csv")
Sys.time()
write.table(myPED, "data/B05_GEMMA/01_PLINK/dpcharMAF20NA10.ped", row.names=FALSE, col.names=FALSE)
Sys.time()

#try removing problematic SNP
myPED[,"39118"]
which( colnames(myPED)=="39118" ) #73399
which( colnames(myPED)=="39118.1" ) #73400
which(myMAP2 == "SNP39118", arr.ind=TRUE)

myPED.ed <- myPED[,-c(73399:73400)]
myMAP.ed <- myMAP2[-c(63),]
write.table(myMAP.ed, "data/B05_GEMMA/01_PLINK/dpcharMAF20NA10_test.map", row.names=FALSE, col.names=FALSE)
Sys.time()
write.table(myPED.ed, "data/B05_GEMMA/01_PLINK/dpcharMAF20NA10_test.ped", row.names=FALSE, col.names=FALSE)
Sys.time()