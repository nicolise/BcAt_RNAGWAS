#Nicole E Soltis
#09/12/18

#--------------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of At reads across ALL GENES 
#from GEMMA_lsm/D_GWAplots 

rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At")
#could redo this with named genes
mydat1 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
mydat10 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top10SNPsample.csv")
mydat100 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top100SNPsample.csv")

##check which file I'm reading in!
mydat <- mydat1

mydat_plot <- mydat[,-c(1)]
names(mydat_plot)[3] <- "chr.snp"

#Make plotting variables for snp
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
mydat_plot$Index.s = NA
lastbase = 0
for (i in unique(mydat_plot$chr.snp)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr.snp==i, ]$Index.s=mydat_plot[mydat_plot$chr.snp==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr.s==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr.s==i, ]$Index.s=mydat_plot[mydat_plot$chr.s==i, ]$ps+lastbase
  }
}

#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$chr.snp)
colScale <- scale_colour_manual(name = "Chromosome",values = myColors)

##double check output name
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/AtCol0_top1SNP_bySNP.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.s, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr.snp),alpha=0.001), stroke=0)+
    labs(list(title=NULL))+
    theme(legend.position="none")+
    #y scale for 1 SNP
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))
)
dev.off()

-log10(min(mydat_plot$p_score))

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=FALSE), "Index.s"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=TRUE), "Index.s"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2


#plot Chr 1 only
mydat_plot_c1 <- mydat_plot[mydat_plot$chr.snp==1,]
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/AtCol0_top1SNP_bySNP_Chr1.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot_c1, aes(x=Index.s, y=(-log10(mydat_plot_c1$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr.snp),alpha=0.001), stroke=0)+
    labs(list(title=NULL))+
    theme(legend.position="none")+
    #y scale for 1 SNP
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Distance (Mb)", breaks=c(0, 1e+6, 2e+6, 3e+6, 4e+6), labels=c(0,1,2,3,4))
)
dev.off()