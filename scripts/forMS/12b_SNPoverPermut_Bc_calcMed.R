#Nicole E Soltis
#07/02/18

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#Bc transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/04_GEMMAout/
#kmat script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A03_GEMMA_kmat.sh
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A04_runGEMMA_kmat_3genos.sh

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/data/")

#get list of phenos with SNP > Thr
mysnplist <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/Bc_phenos_sigSNPovrThr.csv")

#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")
mytime <- Sys.time()
  #each phenotype
j <- "col0"
for (i in mysnplist$pheno){
    #Sys.time()
    my_gemma <- read.table(paste("GEMMA_eachAt_Bc/04_GEMMAout/",j,"/",j,"_MAF20NA10_",i,".assoc.txt", sep=""), header=TRUE)
    #Sys.time()
    my_gemma$pheno <- i
    snpnum5pct <- nrow(my_gemma[my_gemma$p_score < 1.960864e-05,])
    snpnum1pct <- nrow(my_gemma[my_gemma$p_score < 5.022625e-06,])
    pheno <- i
    my_gemma.t <- as.data.frame(cbind(snpnum5pct, snpnum1pct, pheno))
    try(ifelse(i == 8, totnumsnp <- my_gemma.t, totnumsnp <- rbind(totnumsnp, my_gemma.t)))
    }
mytime
Sys.time()


write.csv(totnumsnp, "GEMMA_eachAt_Bc/07_TopSNPs/Bc_TOTAL_numphenosoverThr.csv")

median(totnumsnp$snpnum5pct) #10
median(totnumsnp$snpnum1pct)
blah <- totnumsnp[totnumsnp$snpnum1pct > 0,]
median(blah$snpnum1pct) #13

write.csv(totnumsnp, "07_TopSNPs/Bc_TOTAL_numphenosoverThr.csv")
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
totnumsnp <- read.csv("07_TopSNPs/Bc_TOTAL_numphenosoverThr.csv")

median(totnumsnp$snpnum5pct) #10
median(totnumsnp$snpnum1pct)

#SE = sd / (n^0.5)
mysd <- sd(totnumsnp$snpnum5pct)
myrtn <- nrow(totnumsnp)^0.5
mysd/myrtn

#save as .jpg, 300 wide by 400 tall
hist(log10(totnumsnp$snpnum5pct))
library(ggplot2)
ggplot(totnumsnp, aes(x=log10(snpnum5pct))) + geom_histogram(color="navyblue", fill="royalblue1") + 
  geom_vline(aes(xintercept=median(log10(totnumsnp$snpnum5pct))), color="navyblue", linetype="dashed") + 
  theme_bw()

#for at use darkgreen with palegreen3
#------------------------------------------------------------------------
#for At
rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")

#get list of phenos with SNP > Thr
mysnplistA <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/At_phenos_sigSNPovrThr.csv")

mytime <- Sys.time()
#each phenotype
j <- "col0"
for (i in mysnplistA$pheno){
  #Sys.time()
  my_gemma <- read.table(paste("GEMMA_eachBc_At/04_GEMMAout/col0_round2/",j,"_MAF20NA10_obs_",i,".assoc.txt", sep=""), header=TRUE)
  #Sys.time()
  my_gemma$pheno <- i
  snpnum5pct <- nrow(my_gemma[my_gemma$p_score < 2.901288e-05,])
  snpnum1pct <- nrow(my_gemma[my_gemma$p_score < 6.085872e-06,])
  pheno <- i
  my_gemma.t <- as.data.frame(cbind(snpnum5pct, snpnum1pct, pheno))
  try(ifelse(i == 12, totnumsnp <- my_gemma.t, totnumsnp <- rbind(totnumsnp, my_gemma.t)))
}
mytime
Sys.time()

write.csv(totnumsnp, "GEMMA_eachAt_Bc/07_TopSNPs/At_TOTAL_numphenosoverThr.csv")

setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
totnumsnp <- read.csv("07_TopSNPs/At_TOTAL_numphenosoverThr.csv")

median(totnumsnp$snpnum5pct) #10
median(totnumsnp$snpnum1pct)

#SE = sd / (n^0.5)
mysd <- sd(totnumsnp$snpnum5pct)
myrtn <- nrow(totnumsnp)^0.5
mysd/myrtn

#save as .jpg, 300 wide by 400 tall
hist(log10(totnumsnp$snpnum5pct))
library(ggplot2)
ggplot(totnumsnp, aes(x=log10(snpnum5pct))) + geom_histogram(color="darkgreen", fill="palegreen3") + 
  geom_vline(aes(xintercept=median(log10(totnumsnp$snpnum5pct))), color="darkgreen", linetype="dashed") + 
  theme_bw()

#for at use darkgreen with palegreen3