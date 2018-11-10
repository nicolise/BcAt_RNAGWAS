#Nicole E Soltis
#11/05/18

#-------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of Bc reads across ALL GENES 

#from GEMMA_lsm/D_GWAplots

rm(list=ls())
#start with a small file
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
#for now, don't care about annotating which phenotype was which (gene names for transcripts)
mydat01 <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_RAND2_top1SNPsample.txt")
mydat01 <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_top1SNPsample.csv")

mydat <- mydat01
mydat$chr.snp <- paste(mydat$chr, mydat$ps, sep=".")

#save min p value per permutation
minrow <- mydat[mydat$p_score == min(mydat$p_score),]
#minrowall <- minrow
#do for all 5 permutations
minrowall <- rbind(minrowall, minrow)
minrowall$neg10p <- log10(minrowall$p_score)*-1
#highest = 7.53

#Make plotting variables for snp
mydat_plot <- mydat[order(mydat$chr, mydat$ps),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/5xRand/BcCol0_top1SNP_rand2.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=FALSE), "Index.s"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=TRUE), "Index.s"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
#--------------------------------------------------------------------------------

