#Nicole E Soltis
#07/20/18

#--------------------------------------------------------------------------------------
#manhattan plot

rm(list=ls())
#start with a small file
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At/")
mydat <- read.table("04_GEMMAout/col0_MAF20NA10_obs_1.assoc.txt", sep="\t", header=TRUE)

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

#-------------------------------------------------------------------------------------
#plot by SNP location! 
library(ggplot2)
#create a custom color scale
#myColors <- c("grey20", "grey60")
cols <- c("darkgreen","palegreen3")
myColors <- rep(cols,9)
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/AtCol0_chr1gene1_bySNP.jpg", width=10, height=6, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(chr),alpha=0.001))+
    scale_y_continuous(name="-log10(p)")+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    #geom_hline(yintercept=-log10(0.001), linetype="dashed")+
    expand_limits(y=0)+
    theme(text=element_text(size=16))
)
dev.off()

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=FALSE), "Index.s"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=TRUE), "Index.s"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
