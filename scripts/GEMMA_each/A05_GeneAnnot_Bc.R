#Nicole E Soltis
#07/19/18

#--------------------------------------------------------------------
#this is for any ol Bc genes -- top SNPs genome-wide
#we can also narrow this down for the Bot/ BoA/ Net5/ ABA gene associations

rm(list=ls())
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

#------------------------------------------------------------------------------------------------
##customize for each SNP list
#snp location:
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/06_GEMMAsumm/GeneNames/")
my_files <- list.files(pattern = "\\.csv$")
my_data <- lapply(my_files, read.csv)
names(my_data) <- gsub("\\.csv$", "", my_files)

##will need to come back and redo this for top 100 SNPs!
#------------------------------------------------------------------------------------------------
#read in single SNP list
for (k in 1:length(my_data)){
my.snps <- my_data[[k]]
#calculate gene center
#calculate distance gene center to SNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1
#1:18 but have no sig SNPs on chr 17, 18 so actually 1:16

#associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  my.snp.sub <- my.snps[my.snps$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  for (i in c(1:nrow(my.snp.sub))){
  #for (i in 1:3){
    this.snp <- as.numeric(my.snp.sub[i,"ps"])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    #this prevents erroring out if no SNPs within a certain chromosome
    ifelse(nrow(this.gene) == 0, this.gene <- gtf.sub[1,], this.gene <- this.gene)
    this.line <- cbind(my.snp.sub[i,], this.gene)
    this.line$closest.end <- NA
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, full.snp.genes <- all.genes, full.snp.genes <- rbind(full.snp.genes, all.genes))
  #now only keep genes if nearest end is within +-1kb of SNP (2kb window)
  full.snp.genes$closest.end <- pmin(abs(full.snp.genes$ps - full.snp.genes$V4),abs(full.snp.genes$ps - full.snp.genes$V5)) 
  full.genes.sub <- full.snp.genes[full.snp.genes$closest.end < 1000,]
}
#remove dummy rows when phenotype = NA
full.snp.genes <- full.snp.genes[!is.na(full.snp.genes$pheno),]
assign(paste(names(my_data)[k], "genes", sep="_"),full.snp.genes)
write.csv(full.snp.genes, paste("SNPannot/",names(my_data)[k], "_genes",".csv", sep=""))
}

#full.snp.genes now has all SNPs with gene annotations 

## check file names for threshold
#setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_eachAt_Bc")
#write.csv(full.snp.genes, ".csv")
