#Nicole E Soltis
#020518

#file outputs go to 06_genetoSNPplot_B05.R
#---------------------------------------------------------
library(tidyr)
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/testnames/")
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/outfiles/")
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

#select just all files on Chr1
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/outfiles/")
my.files <- list.files(pattern = "Bcin01g")

#now, keep only top 100 SNPs/gene
pList.100 <- list()
#dummy start, fix this
Chr1.all.top.100 <- as.data.frame(my.file)
Chr1.all.top.100 <- Chr1.all.top.100[1,]
for(i in 1:length(my.files)) {
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  my.file$AbsEst <- abs(my.file[,2])
  my.file$gene <- names(my.file)[2]
  names(my.file)[2] <- "Estimate"
  my.file <- top_n(my.file, 100, AbsEst)
  pList.100[[i]] <- my.file
  Chr1.all.top.100 <- rbind(Chr1.all.top.100, my.file)
}
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
write.csv(Chr1.all.top.100, "Chr1_top100SNPperGene.csv")

#now only keep top 10 per gene
pList.10 <- list()
Chr1.all.top.10 <- my.tops.10
Chr1.all.top.10 <- Chr1.all.top.10[0,]
for (j in 1:length(unique(Chr1.all.top.100$gene))){
#for (j in 1:1){
  my.tops.10 <- Chr1.all.top.100[Chr1.all.top.100$gene==unique(Chr1.all.top.100$gene)[j],]
  my.tops.10 <- top_n(my.tops.10, 10, AbsEst)
  pList.10[[j]] <- my.tops.10
  Chr1.all.top.10 <- rbind(Chr1.all.top.10, my.tops.10)
}

#and top 1 per gene
pList.1 <- list()
Chr1.all.top.1 <- my.tops.1
Chr1.all.top.1 <- Chr1.all.top.1[0,]
for (j in 1:length(unique(Chr1.all.top.10$gene))){
  #for (j in 1:1){
  my.tops.1 <- Chr1.all.top.10[Chr1.all.top.10$gene==unique(Chr1.all.top.10$gene)[j],]
  my.tops.1 <- top_n(my.tops.1, 1, AbsEst)
  pList.1[[j]] <- my.tops.1
  Chr1.all.top.1 <- rbind(Chr1.all.top.1, my.tops.1)
}

#keeps blocks of adjacent SNPs if they have identical estimates -- keep in data for now
Chr1.all.top.1 <- unique(Chr1.all.top.1[,1:4])

#write out both
write.csv(Chr1.all.top.10, "Chr1_top10SNPperGene.csv")
write.csv(Chr1.all.top.1, "Chr1_top1SNPperGene.csv")

#-----------------------------------------------------------------------------------------------------
#now, select just top SNP on all chromosomes (eee!)
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/outfiles/")
my.files <- list.files(pattern = "Bcin")

#now, keep only top 100 SNPs/gene
pList.1 <- list()
#dummy start, fix this
Chr.all.top.1 <- as.data.frame(my.file)
Chr.all.top.1 <- Chr.all.top.100[1,]
for(i in 1:length(my.files)) {
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  my.file$AbsEst <- abs(my.file[,2])
  my.file$gene <- names(my.file)[2]
  names(my.file)[2] <- "Estimate"
  my.file <- top_n(my.file, 1, AbsEst)
  pList.1[[i]] <- my.file
  Chr.all.top.1 <- rbind(Chr.all.top.1, my.file)
}
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
write.csv(Chr.all.top.1, "ChrAll_top1SNPperGene.csv")