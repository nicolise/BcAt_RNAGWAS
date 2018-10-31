#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
#using same genotype input from B05.10 BcAtGWAS bigRR and GEMMA
#removing isolate 01.01.12 prior to GEMMA: GEMMA will drop fewer SNPs
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

setwd("~/Projects/BcGenome/data")
setwd("~/Documents/GitRepos/BcGenome/data")
SNPnames <- read.csv("BO5_97_iso_small/File_key_in_Bo5bamfolder_NES.csv", header=TRUE)
SNPnames <- SNPnames[,c("Isolate","names")]
names(SNPnames)[1]<- "Isolate"
names(SNPs_renamed) <- SNPnames[match(names(SNPs_renamed),SNPnames[,"names"]),"Isolate"] 
SNPs_renamed <- SNPs_renamed[,-c(1)]
names(SNPs_renamed)[1] <- "Chrom"
names(SNPs_renamed)[2] <- "Pos"
SNPs_renamed <- SNPs_renamed[,-c(38)] # remove 1.01.06.1
#reorder columns
SNPs_renamed <- SNPs_renamed[,order(names(SNPs_renamed))]
SNPs_renamed <- SNPs_renamed[,c(66,92,1:65,67:91,93:98)]
#remove 01.01.12
SNPs_renamed <- SNPs_renamed[,c(1:7,9:98)]

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
write.csv(SNPs_renamed, "01_PLINK/OriginalSNPdata.csv")
mySNPs <- SNPs_renamed

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
#MAP2 still has chromosome 1:18

write.table(myMAP2, "01_PLINK/dpbinMAF20NA10.map", row.names=FALSE, col.names=FALSE)

#myMAP2 <- read.table("data/B05_GEMMA_les/E_01_PLINK/dpbinMAF20NA10.map")

#write.csv(mySNPs3, "data/B05_GEMMA_les/D_01_PLINK/dp_binMAF20_10NA.csv")
#write.csv(mySNPs, "data/B05_GEMMA_les/D_01_PLINK/hp_binMAF20_10NA.csv")
Sys.time()
write.table(myPED, "01_PLINK/dpbinMAF20NA10.ped", row.names=FALSE, col.names=FALSE)
write.csv(myPED, "01_PLINK/dpbinMAF20NA10ped.csv")
Sys.time()
