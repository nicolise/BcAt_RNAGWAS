#Nicole E Soltis 
#092217
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------
rm(list=ls())
sink(type=c("output","message"), split=TRUE)
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
SNPs <- read.csv("allreadsGWAS/02_bigRR/lsmeans_allreads_SNPS_MAF20.csv", row.names = 1)
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
Phenos <- read.csv("allreadsGWAS/02_bigRR/lsmeans_allreads_phenos_MAF20.csv", row.names = 1)
#for poisson all values must be positive: additive transformation
#actual phenotypes 2:9268
dat <- Phenos[,c(2:9268)]
min(dat)
dat <- dat + 29

#check which phenos are done
# con <- file("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218.2.csv","r")
# first_line <- readLines(con,n=1)
# close(con)
#see document notes/lsm_bigRR_readme.docx to see which files have which phenotypes
which (colnames(Phenos)=="Bcin03g06050.1") 
which (colnames(dat)=="Bcin14g00260.1")


#--------------------------------------------------------------
#try clearing out memory
#memory.size(max=F)
#gc()
#can also restart R or the machine after a crash

#find current pheno to start from
dat <- dat[,c(6487:9267)]
#make dat column number a multiple of 500
dat[,c(2782:3000)] <- NA

mysplit.dat <- lapply(seq(1,ncol(dat)-500,500), function(u) dat[,u:(u+500)])

# loop from here to run batches of phenotypes in bigRR
for (y in 1:5){
  dat <- mysplit.dat[[y]]
  con <- file(paste("allreadsGWAS/03_bigRRout/lsm_allphenos_MAF20_013118b",y,"log",sep="."))
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
    #testing with impute=T
    Pheno.BLUP.result <- try(bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE, impute = TRUE))
    #can add try here as well
    Pheno.HEM.result <- try(bigRR_update(Pheno.BLUP.result, SNPs))
  
   outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
    colnames(outpt.HEM)[2] <- paste(colnames(dat)[i],"HEM",sep=".")
    #write out to .csv after each phenotype! This saves our progress in case of memory error
    #but is probably slow. going to try a new file for each phenotype
    #write.csv(outpt.HEM, file=paste("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218",y,"csv",sep="."), append=T)
    write.csv(outpt.HEM, file=paste("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_013118b",y,i,"csv",sep="."))
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
}

#-------------------------------------------------------------

#read output back in
partial.out <- read.csv("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218.csv")

library(bigmemory)
system.time(mypartial.out <- read.big.matrix("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218.csv", header = TRUE))

#user  system elapsed 
#208.340   2.464 211.073 
#Warning messages:
#  1: In na.omit(as.integer(firstLineVals)) : NAs introduced by coercion
#2: In na.omit(as.double(firstLineVals)) : NAs introduced by coercion
#3: In read.big.matrix("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218.csv",  : Because type was not specified, we chose double based on the first line of data.

library(data.table)
# system.time(mypartial.out.b <- fread("allreadsGWAS/03_bigRRout/lsm_bigRR_MAF20_012218.csv", header = T, sep = ',')) 
#Read 287826 rows and 920 (of 920) columns from 5.753 GB file in 00:07:20
#user  system elapsed 
#437.240   4.280 441.469 

