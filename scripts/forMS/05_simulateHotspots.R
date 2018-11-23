#Nicole E Soltis

#Simulate random draws as false GEMMA
#check max hotspot
#--------------------------------------------------------------------------------
rm(list=ls())

#check number of SNPs
setwd("~/Projects/BcAt_RNAGWAS")
nrow(read.table("data/GEMMA_eachAt_Bc/01_PLINK/binMAF20NA10.bim"))
#271495 SNPs

#For Bc genes
#check number of genes
nrow(read.table("data/GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt"))
#9268

#make a vector of length 9268. For each element, randomly draw a number from 1:271495. 
#then, draw a histogram of this vector to see the frequency of hits
#trying: replicating this 1000x because why not
mytime <- Sys.time()
for (i in (1:1000)){
  randBc <- as.data.frame(ceiling(runif(9268, min = 0, max = 271495)))
  names(randBc)[1] <- "SNPdraws"
  hist(randBc$SNPdraws)
  randBc$count <- "blah"
  #count # genes with each SNP as top 1 hit
  randBcCt <- aggregate(count ~ SNPdraws, data = randBc, FUN = function(x){NROW(x)})
  randBcCt$Hist <-"blah"
  randBcHist <- aggregate(Hist ~ count, data = randBcCt, FUN = function(x){NROW(x)})
  randBcHist$run <- i 
  ifelse( i == 1, RandRepBc <- randBcHist, RandRepBc <- rbind(RandRepBc, randBcHist))
}
mytime
Sys.time()

table(RandRepBc$count)
#in 15 of 1000 runs, a SNP occurs 4 times in this list
names(RandRepBc) <- c("GenesPerSNP","nSNP","RepRun")
write.csv(RandRepBc, "data/GWA_Sims/Bc_1000x_HotspotOverlap.csv")

RandRepBc <- read.csv("data/GWA_Sims/Bc_1000x_HotspotOverlap.csv")
#violin plot- x is each of these GenesPerSNP levels, y is mean +- variation for nSNP
library(ggplot2)
mybc_plot <- RandRepBc
p <- ggplot(mybc_plot, aes(x=GenesPerSNP, y=nSNP, group=GenesPerSNP))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/SimHotspots/SimHotspots_1000x_BcGenes.jpg", sep=""), width=3, height=4, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="royalblue1", alpha=0.4) +
  labs(x="Number of Genes with SNP hit", y = "Frequency")+ 
  theme_bw()+
  facet_grid(GenesPerSNP ~. , scales='free_y')
dev.off()


#for At genes
#check number of genes
nrow(read.table("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv"))
#23948

mytime <- Sys.time()
for (i in (1:1000)){
  randAt <- as.data.frame(ceiling(runif(23948, min = 0, max = 271495)))
  names(randAt)[1] <- "SNPdraws"
  hist(randAt$SNPdraws)
  randAt$count <- "blah"
  #count # genes with each SNP as top 1 hit
  randAtCt <- aggregate(count ~ SNPdraws, data = randAt, FUN = function(x){NROW(x)})
  randAtCt$Hist <-"blah"
  randAtHist <- aggregate(Hist ~ count, data = randAtCt, FUN = function(x){NROW(x)})
  randAtHist$run <- i 
  ifelse( i == 1, RandRepAt <- randAtHist, RandRepAt <- rbind(RandRepAt, randAtHist))
}
mytime
Sys.time()
table(RandRepBc$GenesPerSNP)

table(RandRepAt$count)
names(RandRepAt) <- c("GenesPerSNP","nSNP","RepRun")
write.csv(RandRepAt, "data/GWA_Sims/At_1000x_HotspotOverlap.csv")
RandRepAt <- read.csv("data/GWA_Sims/At_1000x_HotspotOverlap.csv")
table(RandRepAt$GenesPerSNP)
#violin plot- x is each of these GenesPerSNP levels, y is mean +- variation for nSNP
library(ggplot2)
myat_plot <- RandRepAt
p <- ggplot(myat_plot, aes(x=GenesPerSNP, y=nSNP, group=GenesPerSNP))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/SimHotspots/SimHotspots_1000x_AtGenes.jpg", sep=""), width=3, height=4, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="palegreen3", alpha=0.4) +
  labs(x="Number of Genes with SNP hit", y = "Frequency")+ 
  theme_bw()+
  facet_grid(GenesPerSNP ~. , scales='free_y')
dev.off()