#Nicole E Soltis
#07/02/18

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#Bc transcripts GEMMA
#"~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc/04_GEMMAoutput"
#and 04_GEMMAoutput.tar.gz
#kmat script: data/B05_GEMMA/GEMMA_allphenos_kmatrix.sh
#script: data/B05_GEMMA/runGEMMA_allphenos_kmat.sh
#media/nesoltis/Data/Kliebenstein/Soltis/BcAt_RNAGWAS/B05_GEMMA_Bc/03_fullOutput_nok.tar.gz

#At transcripts GEMMA
#Data/Kliebenstein/Soltis/BcAt_RNAGWAS/B05_GEMMA_At/C05_GEMMAout/all20kAtgene_GEMMA.tar.gz (116 Gb)
#kmat script: Data/Kliebenstein/Soltis/BcAt_RNAGWAS/B05_GEMMA_At/C04_runGEMMA_allAt_kmat.sh
#script: Data/Kliebenstein/Soltis/BcAt_RNAGWAS/B05_GEMMA_At/C04_runGEMMA_allAt_kmat_run.sh

rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc/")


#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")

  #each phenotype
  for (i in 1:9267){
    #actually: 1:9267
    Sys.time()
    my_gemma <- read.table(paste("04_GEMMAoutput/binMAF20NA10_PLINK_",i,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 4 seconds to read 1 phenotype
    #times 9300 = 10.3 hours
    #take top 10 SNP/phenotype
    #also save top 1 SNP/ phenotype
    my_gemma_top100 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top100 <- my_gemma_top100[1:100,]
    my_gemma_top100$pheno <- i
    row.names(my_gemma_top100) <- c(1:100) + (100*(i-1))
    my_gemma_top10 <- my_gemma_top100[1:10,]
    my_gemma_top1 <- my_gemma_top10[1,]
    #this gives an error but it's fine
    try(ifelse( i == 1, write.table(my_gemma_top100, "05_GEMMAsumm/GEMMA_top100SNPsample.txt", sep = ",", col.names = TRUE), write.table(my_gemma_top100, "05_GEMMAsumm/GEMMA_top100SNPsample.txt", sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top10, "05_GEMMAsumm/GEMMA_top10SNPsample.txt", sep = ",", col.names = TRUE), write.table(my_gemma_top10, "05_GEMMAsumm/GEMMA_top10SNPsample.txt", sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top1, "05_GEMMAsumm/GEMMA_top1SNPsample.txt", sep = ",", col.names = TRUE), write.table(my_gemma_top1, "05_GEMMAsumm/GEMMA_top1SNPsample.txt", sep = ",", col.names = FALSE, append = TRUE)))
    Sys.time()
  }

#check range of beta values
mylgSNP <- mytop100[order(-abs(mytop100$beta)),] #0.70 at tops, p ~ 1e-07
mylgSNP <- mytop1[order(-abs(mytop1$beta)),] #0.70 at tops, p ~ 1e-07
head(mylgSNP)

#save just largest effects SNPs - top 26k
hist(abs(mytop100$beta))
mylgsnp <- mytop100[abs(mytop100$beta) > 0.5,]
write.csv(mylgsnp, "data/B05_GEMMA_Bc/05_GEMMAsumm/GEMMA_top100_beta05SNP.csv")

##run this today
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc/")
#now try: z-scale each phenotype and take largest fx SNPs
#each phenotype
for (i in 1:9267){
  #actually: 1:9267
  mystart.time <- Sys.time()
  my_gemma <- read.table(paste("04_GEMMAoutput/binMAF20NA10_PLINK_",i,".assoc.txt", sep=""), header=TRUE)
  Sys.time()
  #takes 4 seconds to read 1 phenotype
  #times 9300 = 10.3 hours
  #try z-scaling 
  my_gemma.z <- my_gemma
  my_gemma.z$beta_z <- scale(my_gemma.z$beta, center = TRUE, scale = TRUE)
  my_gemma.z <- my_gemma.z[abs(my_gemma.z$beta_z) > 4,]
  #this gives an error but it's fine
  try(ifelse( i == 1, write.table(my_gemma.z, "05_GEMMAsumm/GEMMA_topSNPsample_zscale.txt", sep = ",", col.names = TRUE), write.table(my_gemma.z, "05_GEMMAsumm/GEMMA_topSNPsample_zscale.txt", sep = ",", col.names = FALSE, append = TRUE)))
  Sys.time()
}
mystart.time
Sys.time()
#expect approx. 30 SNPs per phenotype (gene)
my_gemma.z.b <- read.table("05_GEMMAsumm/GEMMA_topSNPsample_zscale.txt")
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