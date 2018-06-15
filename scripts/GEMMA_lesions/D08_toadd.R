#Nicole E Soltis
#06/14/18
#---------------------------------------------------------------------------
#fake df with 1 snp 100x
#Does GEMMA do a pairwise search to find SNPs that are 100% identical in their variation and then chose one to use in the model?

#add to previous file:
dropSNPs <- read.csv("B05_GEMMA_les/D_01_PLINK/DroppedSNPs_GEMMA.csv")
justSNPs <- dropSNPs[,4:99]
unique(justSNPs[2,])
justSNPs$count <- apply(justSNPs, 1, function(x)length(unique(x)))
justSNPs$na_count  <- apply(justSNPs, 1, function(x) sum(is.na(x)))
justSNPs$allele_num <- ifelse(justSNPs$na_count > 0, justSNPs$count - 1, justSNPs$count) 
table(justSNPs$allele_num) #cool, all SNPs do in fact have 2 alleles = polymorphic

#pipeline note:
#1. run D01_TABtoPEDnMAP.R
#2. copy plink executable to B05_GEMMA_les/E_01_testSNPs
#3. in command prompt: cd to B05_GEMMA_les/E_01_testSNPs
#4. RUN ./plink --noweb --file run/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE. NEXT STEP is customized by ogphenos/ permutation
#5. run this script (02_prepPhenos.R)
#6. cd to GEMMA_files
#7. copy edited .fam, original .bim, .bed to D_02_randGEMMA/
#8. copy bash script: cp scripts/GEMMA_lesions/norand_GEMMA_kmatrix.sh data/B05_GEMMA_les/
#9. cd to data/B05_GEMMA_les/
#9. calculate k-matrix with: bash norand_GEMMA_kmatrix.sh, mv files to D_03_kmat
#10. run GEMMA: bash norand_GEMMA_kmatrix_run.sh


rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#use same SNP set from B05.10 bigRR
mySNPs <- read.csv("allreads_bigRR/B05.10/01_prepFiles/hp_binMAF20_20NA.csv")

mySNPs.sub <- mySNPs[mySNPs$Chrom==18,]
mySNPs.sub <- mySNPs.sub[mySNPs.sub$Pos > 343603,]

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
write.csv(SNPs_renamed, "B05_GEMMA/01_PLINK/OriginalSNPdata.csv")
mySNPs <- SNPs_renamed
#narrow down to test on only 100 SNPs
mySNPs <- mySNPs[1:100,]

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
#MAP2 still has chromosome 1:18
write.table(myMAP2, "data/B05_GEMMA_les/E_01_testSNPs/dpbinMAF20NA10.map", row.names=FALSE, col.names=FALSE)

myMAP2 <- read.table("data/B05_GEMMA_les/E_01_testSNPs/dpbinMAF20NA10.map")

write.csv(mySNPs3, "data/B05_GEMMA_les/E_01_testSNPs/dp_binMAF20_10NA.csv")
write.csv(mySNPs, "data/B05_GEMMA_les/E_01_testSNPs/hp_binMAF20_10NA.csv")
Sys.time()

#check MAF and missingness threshold
table(myPED[,15]) #this one is fine, except PLINK has trouble finding 0 as an allele before 2... but 0 comes later?! what is up?
myPED[,7:107] <- myPED[,15]
#now every SNP has identical genotype variation to SNP 1

Sys.time()
write.table(myPED, "data/B05_GEMMA_les/E_01_testSNPs/dpbinMAF20NA10.ped", row.names=FALSE, col.names=FALSE)
write.csv(myPED, "data/B05_GEMMA_les/E_01_testSNPs/dpbinMAF20NA10ped.csv")
Sys.time()

#copy fam file over from:
#myFAM_match2 <- read.table("data/B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)

#um, what did I do to fix PLINK being confused about number of alleles?
