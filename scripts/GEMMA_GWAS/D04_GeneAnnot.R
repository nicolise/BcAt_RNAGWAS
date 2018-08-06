#Nicole E Soltis
#07/19/18

#--------------------------------------------------------------------------
#this is for any ol Bc genes -- top SNPs genome-wide
#we can also narrow this down for the Bot/ BoA/ Net5/ ABA gene associations, but need to get those gene IDs first.

##start from here, really

rm(list=ls())
#gtf location: 
#C:/Users/nesol/Documents/Projects/BcGenome/data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf

#snp location:
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")

##from bot/ boa script
#BOA gene annotation
my.boa <- mySNP_boa_named 
(my.boa[2,1]) #first SNP
(my.boa[2,length(my.boa)]) #last SNP
blah <- my.gtf.c1[my.gtf.c1$V4 < 5133,]
max(blah$V4) #none, so start at first SNP:
min(my.gtf.c1$V4)
my.gtf.boa <- my.gtf.c1[my.gtf.c1$V4 > 5428,]
blah <- my.gtf.c1[my.gtf.c1$V5 > 61757,]
min(blah$V5) 
my.gtf.boa <- my.gtf.boa[my.gtf.boa$V5 < 63780,]

my.genes.boa <- my.gtf.boa[,c("V4","V5","V12")]
#get gene start and stop for each
library("plyr")
gene.ends.boa <- ddply(my.genes.boa, c("V12"), summarise,
                       geneMin = min(V4),
                       geneMax = max(V5))

#also need to annotate gene # (GEMMA) to gene names

## start here
num.genes <- my.gtf[unique(my.gtf$V12),]

my.gtf <- my.gtf[,1:14]

#calculate gene center
#calculate distance gene center to SNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1
#1:18 but have no sig SNPs on chr 17, 18

#do need to keep p-scores: 
#for snp-level overlap, have "SUMM" to say how many phenotypes are significantly associated with a particular SNP
#but for gene-level overlap, need to check again in next script (B07) how many phenotypes have a significant association with any SNP in that GENE
plant12.snp <- plant12.snp[,c(2,3,45,58,11,14,17,20,23,26,29,32,35,38,41,44)]

#associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant12.sub <- plant12.snp[plant12.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(plant12.sub))){
    this.snp <- as.numeric(plant12.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant12.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, plant12.genes <- all.genes, plant12.genes <- rbind(plant12.genes, all.genes))
}
#plant12.genes now has all SNPs with gene annotations 

plantHO.snp <- plantHO.snp[,c(2,3,58,45,11,14,17,20,23,26,29,32,35,38,41,44)]
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- plantHO.snp[plantHO.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(plant.sub))){
    this.snp <- as.numeric(plant.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, plantHO.genes <- all.genes, plantHO.genes <- rbind(plantHO.genes, all.genes))
}
#plantHO.genes now has all SNPs with gene annotations

domest.snp <- domest.snp[,c(3,4,11,12,6,8,10)]
#domesticated!
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- domest.snp[domest.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(plant.sub))){
    this.snp <- as.numeric(plant.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, domest.genes <- all.genes, domest.genes <- rbind(domest.genes, all.genes))
}
#domest.genes now has all SNPs with gene annotations

#now only keep genes if nearest end is within 1kb of SNP
domest.genes$closest.end <- pmin(abs(domest.genes$ps - domest.genes$V4),abs(domest.genes$ps - domest.genes$V5)) 
domest.genes.sub <- domest.genes[domest.genes$closest.end < 1000,]

plant12.genes$closest.end <- pmin(abs(plant12.genes$ps - plant12.genes$V4),abs(plant12.genes$ps - plant12.genes$V5)) 
plant12.genes.sub <- plant12.genes[plant12.genes$closest.end < 1000,]

plantHO.genes$closest.end <- pmin(abs(plantHO.genes$ps - plantHO.genes$V4),abs(plantHO.genes$ps - plantHO.genes$V5)) 
plantHO.genes.sub <- plantHO.genes[plantHO.genes$closest.end < 1000,]


## check file names for threshold
write.csv(domest.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_domestgenes_99thr.csv")
write.csv(plantHO.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_plant12HOgenes_99thr.csv")
write.csv(plant12.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_plant12topgenes_99thr.csv")
