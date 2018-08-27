#Nicole E Soltis
#03/05/18

#--------------------------------------------------------------------
#read in GEMMA outputs
#start by focusing on BOT genes?

rm(list=ls())
setwd()
read.csv("data/B05_GEMMA/output")

#read in individual bigRR output files (1 per geno)
#read in all files in folder by pattern matching
my.files <- list.files(pattern = ".csv")

#rename all files
my.filename.key <- as.data.frame(NA)
for(i in 1:length(my.files)) {
  #read only top row
  my.file <- read.csv(my.files[i], nrows=1)
  my.name <- names(my.file)[3]
  file.rename(from=file.path(my.files[i]), to=file.path(paste(my.name,".csv",sep="")))
  my.filename.key[i,1] <- my.files[i]
  my.filename.key[i,2] <- paste(my.name, ".csv", sep="")
}
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