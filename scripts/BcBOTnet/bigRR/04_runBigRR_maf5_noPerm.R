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

## RAN first time WITH POISSON. Lesion size expected to be Gaussian
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
SNPs <- read.csv("03_bigRR/BOT_genos_MAF05.csv", row.names = 1)
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
Phenos <- read.csv("03_bigRR/BOT_phenos_MAF05.csv", row.names = 1)
#try with only 10 phenos first
#total: 2:22
dat <- as.data.frame(Phenos[,c(2:12)])  #INSERT PHENOTYPE COLUMNS HERE

outpt.HEM <- colnames(SNPs)
#thresh.HEM <- list("pos0.95Thresh" = NA, "pos0.975Thresh" = NA, "pos0.99Thresh" = NA, "pos0.999Thresh" = NA, "neg0.95Thresh" = NA, "neg0.975Thresh" = NA, "neg0.99Thresh" = NA, "neg0.999Thresh" = NA)

con <- file("04_bigRRouts/Allphenos_MAF05_092917.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

#Calculate HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1) #good to here
  
  #added try here
  #testing with impute=T
  Pheno.BLUP.result <- try(bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE, impute=TRUE)) #why is this failing for Col0.AT3G26830 and Col0.AT2G30770
  #can add try here as well
  Pheno.HEM.result <- try(bigRR_update(Pheno.BLUP.result, SNPs))
  
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  
  #Permute Thresholds for Phenos - this is what takes forever
  #perm.u.HEM <- vector()
  #for(p in 1:1000) {
  # if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
  # try(temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE))
  # try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE, impute=TRUE),silent = TRUE)
  # try(temp.HEM <- bigRR_update(temp.BLUP, SNPs)) #REF change- was bigRR_update(Pheno.BLUP.result...
  # perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
  
  #}
  #write.csv(perm.u.HEM, paste("data/04_bigRRouts/PermEffects_",colnames(dat)[i],".csv",sep=""))
  #thresh.HEM$"pos0.95Thresh"[i] <- quantile(perm.u.HEM,0.95)
  #thresh.HEM$"pos0.975Thresh"[i] <- quantile(perm.u.HEM,0.975)
  #thresh.HEM$"pos0.99Thresh"[i] <- quantile(perm.u.HEM,0.99)
  #thresh.HEM$"pos0.999Thresh"[i] <- quantile(perm.u.HEM,0.999)
  #thresh.HEM$"neg0.95Thresh"[i] <- quantile(perm.u.HEM,0.05)
  #thresh.HEM$"neg0.975Thresh"[i] <- quantile(perm.u.HEM,0.025)
  #thresh.HEM$"neg0.99Thresh"[i] <- quantile(perm.u.HEM,0.01)
  #thresh.HEM$"neg0.999Thresh"[i] <- quantile(perm.u.HEM,0.001)
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
  
}

# Restore output to console
sink() 
sink(type="message")

#Give column names to the thresholds from the HEM list
#for(j in 1:length(thresh.HEM)) {
#  names(thresh.HEM[[j]]) <- colnames(dat)
#}

#RF-give row names to thresh.HEM and thresh.BLUP so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
#thresh.HEM$"pos0.95Thresh" <- c("pos 0.95 Thresh", thresh.HEM$"pos0.95Thresh")
#thresh.HEM$"pos0.975Thresh" <- c("pos 0.975 Thresh", thresh.HEM$"pos0.975Thresh")
#thresh.HEM$"pos0.99Thresh" <- c("pos 0.99 Thresh", thresh.HEM$"pos0.99Thresh")
#thresh.HEM$"pos0.999Thresh" <- c("pos 0.999 Thresh", thresh.HEM$"pos0.999Thresh")
#thresh.HEM$"neg0.95Thresh" <- c("neg 0.95 Thresh", thresh.HEM$"neg0.95Thresh")
#thresh.HEM$"neg0.975Thresh" <- c("neg 0.975 Thresh", thresh.HEM$"neg0.975Thresh")
#thresh.HEM$"neg0.99Thresh" <- c("neg 0.99 Thresh", thresh.HEM$"neg0.99Thresh")
#thresh.HEM$"neg0.999Thresh" <- c("neg 0.999 Thresh", thresh.HEM$"neg0.999Thresh")

#Write results to output
#write.csv(rbind(thresh.HEM$"pos0.95Thresh",thresh.HEM$"pos0.975Thresh",thresh.HEM$"pos0.99Thresh",thresh.HEM$"pos0.999Thresh",thresh.HEM$"neg0.95Thresh",thresh.HEM$"neg0.975Thresh",thresh.HEM$"neg0.99Thresh",thresh.HEM$"neg0.999Thresh",outpt.HEM),"04_bigRRouts/bigRR_MAF10_092917.HEM.csv")
write.csv(outpt.HEM,"04_bigRRouts/bigRR_MAF05_092917.HEM.csv")