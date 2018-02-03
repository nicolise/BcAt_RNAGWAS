#Nicole E Soltis 
#092217
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")

#########################
# This makes the bigRR_update run through the GPU
# You need to do this first to mask the native 'bigRR_update' in the bigRR package
# one alternative to family = gaussian(link = identity) is family = poisson(link = log)

## Lesion size is Gaussian (Wei Zhang 2018)
## RNAseq should be Poisson
bigRR_update <- function (obj, Z, family = poisson(link = identity), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  #if bigRR is having trouble with missing values (NAs) can add option impute=TRUE
  #X is the genotype (MyX)
  #y is the phenotype (dat)
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE, impute = TRUE )
}
########################
#NOTE1 FROM RACHEL:  we need bigRR1.3-9 to get GPU option
# must download R-Forge version bigRR1.3-9tar.gz and manually install
# https://r-forge.r-project.org/R/?group_id=1301
# install.packages("bigRR", repos="http://R-Forge.R-project.org")

#NOTE2 FROM RACHEL: need package 'gputools' but CRAN version fails to install
# must first install Nvidia's CUDA toolkit -current version is 7.5
# installed from developer.nvidia.com/cuda-downloads

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("allreadsGWAS/02_bigRR/col0_allreads_SNPS_MAF20.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Con.Pos <- do.call(paste, c(SNPs[c("Chromosome","Contig","Pos")], sep="."))

rownames(SNPs) <- SNPs[,"Chr.Con.Pos"] #set the new column of chrom.base as rownames 
any(duplicated(SNPs$Chr.Con.Pos))#check that none are duplicated
SNPs <- SNPs[,4:98] #take out first three cols (X.CHROM, POS, REF) and new last column
ogSNPs <- SNPs

#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

#read in phenotype data
Phenos <- read.csv("allreadsGWAS/02_bigRR/col0_allreads_phenos_MAF20.csv", row.names = 1)
#for poisson all values must be positive
#actual phenotypes 3:9269
dat <- Phenos[,c(3:9269)]
min(dat)
dat <- dat + 49

#111818 ran to "Bcin11g01190.1"
#which( colnames(dat)=="Bcin11g01190.1" )
#111919
#dat <- dat[,c(586:9267)]

outpt.HEM <- colnames(SNPs)

con <- file("allreadsGWAS/03_bigRRout/col0_allphenos_MAF20_011918.log")
sink(con, append=TRUE, type="message")

print(Sys.time())
#Calculate HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1) #good to here
  
  #added try here
  #testing with impute=T
  Pheno.BLUP.result <- try(bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE, impute=TRUE))
  #can add try here as well
  Pheno.HEM.result <- try(bigRR_update(Pheno.BLUP.result, SNPs))
  
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
  #write out to .csv after each phenotype! This saves our progress in case of memory error
  write.csv(outpt.HEM, file="allreadsGWAS/03_bigRRout/col0_bigRR_MAF20_011918.csv", append=T)
}

# Restore output to console
#sink() 
sink(type="message")
print(Sys.time())