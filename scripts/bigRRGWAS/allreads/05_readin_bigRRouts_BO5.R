#Nicole E Soltis
#020518
#---------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
#read in individual bigRR output files (1 per geno)
#read in all files in folder by pattern matching

#testing: just files on Chr1

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
my.full.data<-my.full.data[,order(colnames(my.full.data),decreasing=TRUE)]
names(my.full.data)

#get gff3 file info
library("ape")
setwd("~/Projects/BcGenome/data/ensembl/BO5.10")
#can extract the files from .gz using 7-zip
#then read using R!
my.gff <- read.gff("extracted/Botrytis_cinerea.ASM83294v1.38.chromosome.1.gff3/Botrytis_cinerea.ASM83294v1.38.chromosome.1.gff3", na.strings = c(".", "?"))
my.gene.list <- as.data.frame(names(my.full.data)[2:31])
names(my.gene.list)[1] <- "mygenes"
my.gene.list$mygenes <- as.character(my.gene.list$mygenes)
#some lazy regex
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$mygenes, ".1.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".2.HEM"), "[", 1)
#now loop over to keep things
#1:length(my.gene.list$justgenes)
my.gff.genes <- my.gff[1,]
my.gff.genes$transcript <- NA
my.gff.genes <- my.gff.genes[-c(1),]
for (y in (1:length(my.gene.list$justgenes))){
  my.current.gene <- my.gff[grep(my.gene.list[y,2], my.gff$attributes), ]
  my.current.gene$transcript <- my.gene.list[y,2]
  my.gff.genes <- rbind(my.gff.genes, my.current.gene)
}

#for a first go, just keep the "gene" rows
my.transcript.locs <- my.gff.genes[my.gff.genes$type=="gene", ]

#now: plot for each transcript
