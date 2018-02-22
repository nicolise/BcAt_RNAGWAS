#Nicole E Soltis 
#092217
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------
rm(list=ls())
sink(type=c("output","message"), split=TRUE)
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")

#########################
# This makes the bigRR_update run through the GPU
# You need to do this first to mask the native 'bigRR_update' in the bigRR package
# one alternative to family = gaussian(link = identity) is family = poisson(link = log)

## Lesion size is Gaussian (Wei Zhang 2018)
bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
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
SNPs <- read.csv("data/testMethods/BO5tomato/02_bigRR/tomato_zscale_SNPS_MAF20.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Pos <- do.call(paste, c(SNPs[c("Chrom","Pos")], sep="."))

rownames(SNPs) <- SNPs[,"Chr.Pos"] #set the new column of chrom.base as rownames 
any(duplicated(SNPs$Chr.Pos))#check that none are duplicated
SNPs <- SNPs[,3:95] #take out first two cols (X.CHROM, POS) and new last column
ogSNPs <- SNPs

#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

#read in phenotype data
Phenos <- read.csv("data/testMethods/BO5tomato/02_bigRR/tomato_zscale_phenos_MAF20.csv", row.names = 1)
#for poisson all values must be positive: additive transformation
#actual phenotypes 2:9268
#check which phenos are done
#see document notes/lsm_bigRR_readme.docx to see which files have which phenotypes
#which (colnames(Phenos)=="Bcin02g08990.1") 
dat <- Phenos[,c(2:16)]

#--------------------------------------------------------------
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
# loop from here to run batches of phenotypes in bigRR
  con <- file(paste("data/testMethods/BO5tomato/03_bigRRout/B05.10_tomato_zscale_bigRR.log"))
  #type = "message"
  sink(con, append=TRUE, split=TRUE)
  time1 <- print(Sys.time())
  #Calculate HEMs for all phenotypes
  for(i in 1:dim(dat)[2]) { #i will be each isolate
    print(colnames(dat)[i])
    outpt.HEM <- NULL
    outpt.HEM <- colnames(SNPs)
    MyX <- matrix(1, dim(dat)[1], 1) #good to here
    
    #added try here
    #use impute=T for missing data
    Pheno.BLUP.result <- try(bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE, impute = TRUE))
    #can add try here as well
    Pheno.HEM.result <- try(bigRR_update(Pheno.BLUP.result, SNPs))
    
    outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
    colnames(outpt.HEM)[2] <- paste(colnames(dat)[i],"HEM",sep=".")
    #write out to .csv after each phenotype! This saves our progress in case of memory error
    #but is probably slow. going to try a new file for each phenotype
    #write.csv(outpt.HEM, file=paste("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218",y,"csv",sep="."), append=T)
    write.csv(outpt.HEM, file=paste("data/testMethods/BO5tomato/03_bigRRout/B05.10_tomato_zscale_bigRR",i,"csv",sep="."))
    print(Sys.time())
    #run garbage collection just in case to free up space
    gc()
  }
  
  # Restore output to console
  print(time1)
  print(Sys.time())
  sink(type="message")
  sink(file=NULL) 
  sink()