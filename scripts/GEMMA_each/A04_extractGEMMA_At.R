#Nicole E Soltis
#08/

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#At transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/output
#kmat script: 
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/A04_runGEMMA_kmat_3genos.sh

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At")

#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")
    
  #each phenotype
for (j in c("col0","coi1","npr1")){
 for (i in 1:23947){
    #actually: 1:23947 for At
    Sys.time()
    my_gemma <- read.table(paste("output/",j,"/",j,"_MAF20NA10_",i,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 4 seconds to read 1 phenotype
    #times 24000 = 27 hours
    #take top 10 SNP/phenotype
    #also save top 1 SNP/ phenotype
    my_gemma$pheno <- i
    my_gemma_top100 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top100 <- my_gemma_top100[1:100,]
    row.names(my_gemma_top100) <- c(1:100) + (100*(i-1))
    my_gemma_top10 <- my_gemma_top100[1:10,]
    my_gemma_top1 <- my_gemma_top10[1,]
    #and z scaling
    my_gemma.z <- my_gemma
    my_gemma.z$beta_z <- scale(my_gemma.z$beta, center = TRUE, scale = TRUE)
    my_gemma.z <- my_gemma.z[abs(my_gemma.z$beta_z) > 4,]
    mylgsnp <- my_gemma[abs(my_gemma$beta) > 0.5,]
    #this gives an error but it's fine
    try(ifelse( i == 1, write.table(my_gemma_top100, paste("05_GEMMAsumm/",j,"_GEMMA_top100SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top100, paste("05_GEMMAsumm/",j,"_GEMMA_top100SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top10, paste("05_GEMMAsumm/",j,"_GEMMA_top10SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top10, paste("05_GEMMAsumm/",j,"_GEMMA_top10SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top1, paste("05_GEMMAsumm/",j,"_GEMMA_top1SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top1, paste("05_GEMMAsumm/",j,"_GEMMA_top1SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma.z, paste("05_GEMMAsumm/",j,"_GEMMA_topSNPsample_zscale.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma.z, paste("05_GEMMAsumm/",j,"_GEMMA_topSNPsample_zscale.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(mylgsnp, paste("05_GEMMAsumm/",j,"_GEMMA_top100_beta05SNP.txt", sep=""), sep = ",", col.names = TRUE), write.table(mylgsnp, paste("05_GEMMAsumm/",j,"_GEMMA_top100_beta05SNP.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    Sys.time()
  }
}
top1_SNP <- read.table("C06_GEMMAsumm/GEMMA_top1SNPsample.txt", sep = ",")
#------------------------------------------------------------------------------------
#more stuff here - extra code
my_gemma_top10b <- read.table("05_GEMMAsumm/GEMMA_top10SNPsample.txt", sep=",")

#now cp D_07_randSUMM to ~/Documents/GitRepos/BcSolGWAS/
#and make sure nesoltis user account has rwx permission for the new files

#double check that things worked!
blah <- read.csv(paste("D_06_randOUT/quantiles/rand1k_",i,"/pheno",j,".csv", sep=""))
blah <- read.csv("D_07_randSUMM/GEMMA_1krand_SNPsample.csv")


#-----------------------------------------------------
#rename all files
my.filename.key <- as.data.frame(NA)
#for(i in 1:length(my.files)) {
i <- 1
  #read only top row
  my.file <- read.table(my.files[i], nrows=1)
  my.name <- names(my.file)[3]
  file.rename(from=file.path(my.files[i]), to=file.path(paste(my.name,".csv",sep="")))
  my.filename.key[i,1] <- my.files[i]
  my.filename.key[i,2] <- paste(my.name, ".csv", sep="")
#}
names(my.filename.key)[1]<- "outputFile"
names(my.filename.key)[2]<- "TranscriptFile"
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/")
write.csv(my.filename.key, "Key_filenames.csv")

#read in a subset of files, extract top 100 SNPS as in bigRR ==> compare
#now, select just top SNP on all chromosomes (eee!)
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/outfiles/")
my.files <- list.files(pattern = "Bcin")

#now, keep only top 1 SNPs/gene
pList.1 <- list()
#dummy start, fix this
Sys.time()
for(i in 2:length(my.files)) {
  #for (i in 1:1){
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  my.file$AbsEst <- abs(my.file[,2])
  my.file$gene <- names(my.file)[2]
  names(my.file)[2] <- "Estimate"
  my.file <- top_n(my.file, 1, AbsEst)
  pList.1[[i]] <- my.file
  #Chr.all.top.1 <- my.file
  Chr.all.top.1 <- rbind(Chr.all.top.1, my.file)
}
Sys.time()
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
write.csv(Chr.all.top.1, "ChrAll_top1SNPperGene.csv")
Sys.time()

#do this again for top 10! fun!
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/outfiles/")
my.files <- list.files(pattern = "Bcin")

#now, keep only top 1 SNPs/gene
pList.10 <- list()
#dummy start, fix this
Sys.time()
for(i in 2:length(my.files)) {
  #for (i in 1:1){
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  my.file$AbsEst <- abs(my.file[,2])
  my.file$gene <- names(my.file)[2]
  names(my.file)[2] <- "Estimate"
  my.file <- top_n(my.file, 10, AbsEst)
  pList.10[[i]] <- my.file
  #Chr.all.top.10 <- my.file
  Chr.all.top.10 <- rbind(Chr.all.top.10, my.file)
}
Sys.time()
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
write.csv(Chr.all.top.10, "ChrAll_top10SNPperGene.csv")
Sys.time()

#---------------------------------------------------------------------------------------------------------------
#if want top 1% SNPs/ pheno
#ntop1 <- round(nrow(my_gemma)*0.01)
#my_gemma_top1 <- my_gemma[order(my_gemma$p_score),]
#my_gemma_top1 <- my_gemma_top1[1:ntop1,]

#if want quantiles
#need to extract
# myquantsdf <- as.data.frame(NULL)
# for (myquant in seq(0.01,1,0.01)){
#   quantout <- quantile(my_gemma$p_score, myquant)
#   #print(quantout)
#   myrow <- myquant*100
#   myquantsdf[myrow,1] <- myrow
#   names(myquantsdf)[1] <- "Quantile"
#   myquantsdf[myrow,2] <- quantout
#   names(myquantsdf)[2] <- "p_score"
#   myquantsdf[myrow,3] <- i
#   names(myquantsdf)[3] <- "phenotype"