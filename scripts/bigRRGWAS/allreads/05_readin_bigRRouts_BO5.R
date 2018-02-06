#Nicole E Soltis
#020518
#---------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5/03_bigRRout")
#read in individual bigRR output files (1 per geno)
#read in all files in folder by pattern matching
my.files <- list.files(pattern = ".csv")
pList <- list()
for(i in 1:length(my.files)) {
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  pList[[i]] <- my.file
}

#merge all phenotypes into one file -- may only make sense for small sets
my.full.data <- pList[[1]]
for(j in 2:length(pList)){
  my.full.data <- cbind(my.full.data, pList[[j]][,2])
  names(my.full.data)[j+1] <- names(pList[[j]])[2]
}

#so: do top SNP locations correlate with gene location? (Cis fx)
#need: key for gene locations from Vivian
#first step: plot top 100 SNPs per transcript
my.full.data.b<-my.full.data[,order(colnames(my.full.data),decreasing=TRUE)]
names(my.full.data)
