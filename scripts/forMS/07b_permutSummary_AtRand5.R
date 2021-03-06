#Nicole E Soltis
#11/05/18

#-----------------------------------------------------------------------------
#At max overlap is 80, on Chromosome 9_974817 -- set to 150 for sig hotspots?
#and 20 for Bc?
#at these levels, we get ~4 hotspots for At and 4 for Bc

#more permissive threshold would be 5 for Bc, 20 for At, then remove all sites with a random hotspot over these thresholds

#summarize across 5 permutations, without thresholding
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
for(i in c(1:5)){
myranddat <- read.csv(paste0("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=''))

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
jpeg(paste0("plots/Manhattans/5xRand/BcCol0_At_RandHots_run",i,"_top1SNP.jpg", sep=""), width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=pheno))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Number of Genes", breaks=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85))+
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
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(fullsumm, "data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/HotspotSumm_5xRand.csv")
write.csv(fullpeaks, "data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
write.csv(full100peaks, "data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_100top.csv")

full100peaks <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_100top.csv")
fullpeaks <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
hipeaks <- fullpeaks[fullpeaks$numGenes > 10,]
write.csv(hipeaks, "data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/PeaksOver10.csv")

#-------------------------------------------------------------------------------------------
#look at permut summary data
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
fullpeaks <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/TopHotspots_3genepeaks.csv")
hipeaks <- fullpeaks[fullpeaks$numGenes > 10,]
fullsumm <- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/HotspotSumm_5xRand.csv")

i  <- 4
myranddat <- read.csv(paste0("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=''))
#summary info of p-values
quantile(myranddat$p_score, c(0.00, 0.01, 0.05, 0.6))

#          0%           1%           5% 
#3.461793e-08 6.564719e-06 3.014453e-05
#4.923707e-08 5.760712e-06 2.742912e-05 
#4.211615e-08 6.002506e-06 2.886701e-05
#5.873342e-08 5.530269e-06 2.790771e-05 
#2.926053e-08 6.571152e-06 3.071604e-05
#5% level
mean(c(3.014453e-05, 2.742912e-05, 2.886701e-05, 2.790771e-05, 3.071604e-05))
#1% level
mean(c(6.564719e-06,5.760712e-06,6.002506e-06,5.530269e-06,6.571152e-06))
#real data
#setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS")
#mydat100 <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top100SNPsample_genes.csv")
mydat <- read.csv("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes.csv")
quantile(mydat$p_score, c(0.00, 0.01, 0.05, 0.6))

#on linux
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
mydat100 <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/col0_round1/col0_GEMMA_top100SNPsample.txt")

mydat100 <- read.csv("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top100SNPsample_genes.csv")
hi100 <- mydat100[mydat100$p_score < 2.901288e-05,]
hisumm <- as.data.frame(table(hi100$pheno))
sigsumm <- hisumm[hisumm$Freq > 0,]
5213/28373

topsumm <- hisumm[hisumm$Freq > 99,]
names(topsumm)[1] <- "pheno"
names(sigsumm)[1] <- "pheno"
write.csv(topsumm, "GEMMA_eachAt_Bc/07_TopSNPs/At_phenos_manySNPovrThr.csv")
write.csv(sigsumm, "GEMMA_eachAt_Bc/07_TopSNPs/At_phenos_sigSNPovrThr.csv")
