#Nicole E Soltis
#11/05/18

#-----------------------------------------------------------------------------
#summarize across 5 permutations, without thresholding
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
for(i in c(1:5)){
myranddat <- read.csv(paste0("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=''))

#hotspots: count # genes with each SNP as top 1 hit
#chr.ps is bad: R interprets this as a number and rounds it!
myranddat$chr_ps <- paste(myranddat$chr, myranddat$ps, sep="_")
mydat_hots <- aggregate(pheno ~ chr_ps, data = myranddat, FUN = function(x){NROW(x)})
mydat_labels <- myranddat[,c("chr","ps","chr_ps")]
mydat_plot <- merge(mydat_hots, mydat_labels, by="chr_ps", all=FALSE)
mydat_plot <- unique(mydat_plot)

#plot top rand hotspots!
#Make plotting variables for snp
mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (j in unique(mydat_plot$chr)) {
  print(j)
  if (j==1) {
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==j-1)$ps, 1)
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps+lastbase
  }
}

#plot by SNP location!
library(ggplot2)
#create a custom color scale
myColors <- (rep(c("darkred", "indianred1"), 9))
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste0("plots/Manhattans/5xRand/BcCol0_RandHots_run",i,"_top1SNP.jpg", sep=""), width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=pheno))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Number of Genes", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))
)
dev.off()

mysumm <- as.data.frame(table(mydat_hots$pheno))
names(mysumm)<- c("GeneCount", paste("Freq_run",i,sep=""))
ifelse(i == 1, fullsumm <- mysumm, fullsumm <- merge(fullsumm, mysumm, by="GeneCount", all=TRUE))

#check hotspots peak locations -- take all loci with #genes > 2
hotpeaks <- mydat_plot[mydat_plot$pheno > 2,]
hotpeaks$run <- i
ifelse(i == 1, fullpeaks <- hotpeaks, fullpeaks <- rbind(fullpeaks, hotpeaks))

#also take top 100 hotspots per permutation (arbitrary chunk in whatever the lowest Gene # is)
hot100peaks <- mydat_plot[order(-1*mydat_plot$pheno),]
hot100peaks$run <- i
hot100peaks <- hot100peaks[1:100,]
ifelse(i == 1, full100peaks <- hot100peaks, full100peaks <- rbind(full100peaks, hot100peaks))
}
names(fullpeaks)[2] <- "numGenes"
names(full100peaks)[2] <- "numGenes"
write.csv(fullsumm, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/HotspotSumm_5xRand.csv")
write.csv(fullpeaks, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
write.csv(full100peaks, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/TopHotspots_100top.csv")

fullpeaks <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
hipeaks <- fullpeaks[fullpeaks$numGenes > 5,]
write.csv(hipeaks, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/PeaksOver5.csv")

#---------------------------------------------------------------------------
#next, look at deeper hotspots: plot top 10 SNPs per gene per permuation
#summarize across 5 permutations, without thresholding
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
for(i in c(1:5)){
  myranddat <- read.csv(paste0("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top10SNPsample.txt", sep=''))
  
  #hotspots: count # genes with each SNP as top 1 hit
  #chr.ps is bad: R interprets this as a number and rounds it!
  myranddat$chr_ps <- paste(myranddat$chr, myranddat$ps, sep="_")
  mydat_hots <- aggregate(pheno ~ chr_ps, data = myranddat, FUN = function(x){NROW(x)})
  mydat_labels <- myranddat[,c("chr","ps","chr_ps")]
  mydat_plot <- merge(mydat_hots, mydat_labels, by="chr_ps", all=FALSE)
  mydat_plot <- unique(mydat_plot)
  
  #plot top rand hotspots!
  #Make plotting variables for snp
  mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
  mydat_plot$Index = NA
  lastbase = 0
  #Redo the positions to make them sequential		-- accurate position indexing
  for (j in unique(mydat_plot$chr)) {
    print(j)
    if (j==1) {
      mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps
    }	else {
      lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==j-1)$ps, 1)
      mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps+lastbase
    }
  }
  # #plot by SNP location!
  # library(ggplot2)
  # #create a custom color scale
  # myColors <- (rep(c("darkred", "indianred1"), 9))
  # names(myColors) <- levels(mydat_plot$chr)
  # colScale <- scale_colour_manual(name = "Chrom",values = myColors)
  # 
  # setwd("~/Projects/BcAt_RNAGWAS")
  # jpeg(paste0("plots/Manhattans/5xRand/BcCol0_RandHots_run",i,"_top10SNP.jpg", sep=""), width=8, height=5, units='in', res=600)
  # print(
  #   ggplot(mydat_plot, aes(x=Index, y=pheno))+
  #     theme_bw()+
  #     colScale+
  #     #used stroke = 0 for top 10, not top 1
  #     #, stroke=0
  #     geom_point(aes(color = factor(chr),alpha=0.001))+
  #     labs(list( title=NULL))+
  #     theme(legend.position="none")+
  #     scale_y_continuous(name="Number of Genes", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))+
  #     scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))
  #  )
  # dev.off()
  # 
  # mysumm <- as.data.frame(table(mydat_hots$pheno))
  # names(mysumm)<- c("GeneCount", paste("Freq_run",i,sep=""))
  # ifelse(i == 1, fullsumm <- mysumm, fullsumm <- merge(fullsumm, mysumm, by="GeneCount", all=TRUE))

  #check hotspots peak locations -- take all loci with #genes > 2
  hotpeaks <- mydat_plot[mydat_plot$pheno > 3,]
  hotpeaks$run <- i
  ifelse(i == 1, fullpeaks <- hotpeaks, fullpeaks <- rbind(fullpeaks, hotpeaks))

   # #also take top 100 hotspots per permutation (arbitrary chunk in whatever the lowest Gene # is)
   # hot100peaks <- mydat_plot[order(-1*mydat_plot$pheno),]
   # hot100peaks$run <- i
   # hot100peaks <- hot100peaks[1:100,]
   # ifelse(i == 1, full100peaks <- hot100peaks, full100peaks <- rbind(full100peaks, hot100peaks))
}
names(fullpeaks)[2] <- "numGenes"
names(full100peaks)[2] <- "numGenes"
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(fullsumm, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/SNP10_HotspotSumm_5xRand.csv")
write.csv(fullpeaks, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/SNP10_TopHotspots_3genepeaks.csv")
write.csv(full100peaks, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/SNP10_TopHotspots_100top.csv")

#-------------------------------------------------------------------------------
#for each permutation hotspot, plot # genes at 1-SNP level vs. # genes at 10-SNP level
#try for fullpeaks first
fullpeak1 <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
fullpeak10 <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/SNP10_TopHotspots_3genepeaks.csv")

head(fullpeak10)
fullpeak1 <- fullpeak1[,c(2:5,7)]
names(fullpeak1) <- c("chr_ps","numGenes.1","chr","ps","run.1")
fullpeak10 <- fullpeak10[,c(2:5,7)]
names(fullpeak10) <- c("chr_ps","numGenes.10","chr.10","ps.10","run.10")
peakmerge <- merge(fullpeak1, fullpeak10, by="chr_ps", all=FALSE)

#Make plotting variables for snp
mydat_plot <- peakmerge
mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (j in unique(mydat_plot$chr)) {
  print(j)
  if (j==1) {
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==j-1)$ps, 1)
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$ps+lastbase
  }
}

#plotting these multiple peaks per SNP threshold is nothing. Can also try plotting mean # genes per SNP in this overlap list-- takes into account how many of the permutations found it as a hotspot

peaksumm_plot <- aggregate(cbind(mydat_plot$numGenes.1, mydat_plot$numGenes.10), by=list(Chr_ps=mydat_plot$chr_ps), FUN=sum)
names(peaksumm_plot) <- c("Chr_ps","sumGene.1","sumGene.10")
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(peaksumm_plot, "data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/RandHotspots_Correlation.csv")

library(ggplot2)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste0("plots/Manhattans/5xRand/BcCol0_RandHotspots_Corr.jpg", sep=""), width=8, height=5, units='in', res=600)
print(
ggplot(peaksumm_plot, aes(x=sumGene.1, y=sumGene.10))+
  theme_bw()+
  geom_point(aes(color = factor(Chr_ps),alpha=0.001))+
  labs(list( title=NULL))+
  theme(legend.position="none")+
  scale_y_continuous(name="Number of Genes Associated at Top 10 SNP Level")+
  scale_x_continuous(name="Number of Genes Associated at Top 1 SNP Level")
)
dev.off()