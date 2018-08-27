#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
#using same genotype input from B05.10 BcAtGWAS bigRR and GEMMA
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#use same SNP set from B05.10 bigRR
mySNPs <- read.csv("allreads_bigRR/B05.10/01_prepFiles/hp_binMAF20_20NA.csv")
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
mySNPs <- mySNPs[,-c(38,53)] #remove 1.01.06.1 (duplicated) and 01.01.12 (no phenotypes)

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
mySNPs2 <- mySNPs[,-c(1:2)]

#for PED, NA must be replaced with 0 for genotypes, else NA will be read as an allele
#so first, set all genotypes = 0 to =2
mySNPs2[mySNPs2==0] <- 2
mySNPs2[is.na(mySNPs2)] <- 0

#turn all SNPs to "diploid"
#haha, it takes 4 days to do this as a "for" loop (for each row, rbind twice)
#because is.na <-0 before this step, there should be NO heterozygous SNP calls
#this is super fast:
mySNPs3 <- mySNPs2[rep(1:nrow(mySNPs2),each=2),] 

#transpose and format for PED
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
myMAP <- mySNPs[,c("Chrom","Pos")]
myMAP2 <- myMAP
myMAP2$SNPID <- paste("SNP",myMAP2$Pos, sep="")
myMAP2$SNPcM <- 0
myMAP2 <- myMAP2[,c(1,3,4,2)]
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
write.table(myMAP2, "data/B05_GEMMA/C02_PLINK/dpbinMAF20NA10.map", row.names=FALSE, col.names=FALSE)
#add a column of "SNP identifiers" in excel?

#write.csv(mySNPs3, "data/GEMMA_files/C02_PLINK/dp_charMAF20_10NA.csv")
#write.csv(mySNPs, "data/GEMMA_files/C02_PLINK/hp_charMAF20_10NA.csv")
Sys.time()
write.table(myPED, "data/B05_GEMMA/C02_PLINK/dpbinMAF20NA10.ped", row.names=FALSE, col.names=FALSE)
Sys.time()

#check fam files
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
famCSV <- read.csv("B05_GEMMA_At/C03_runGEMMA/binMAF20NA10_fam.csv")
myfam <- read.table("B05_GEMMA_At/C03_runGEMMA/binMAF20NA10.fam")
