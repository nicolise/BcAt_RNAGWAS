#Nicole E Soltis
#11/05/18

#-----------------------------------------------------------------------------
#will later add thresholding info here. For now, just looking at gene overlap on SNPs

rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At")
mydat01 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
#Make plotting variables for snp
mydat <- mydat01
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

#summarize within each SNP - # of transcript hits
mydat_summ <- mydat_plot[,c("chr","ps","p_score","Gene","Index")]
mydat_summ_ngene <- aggregate(Gene ~ Index, data = mydat_summ, FUN = function(x){NROW(x)})
#now add SNP data back on, matching by Index.s
mydat_summ_ngene <- merge(mydat_summ_ngene, mydat_summ[,c("chr","ps","Index")], by="Index")
#remove duplicate rows
mydat_summ_ngene <- unique(mydat_summ_ngene)
mydat_plot <- mydat_summ_ngene

library(ggplot2)
#create a custom color scale
myColors <- rep(c("darkgreen", "palegreen3"), 9)
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chr",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
## check name depending on df
jpeg("plots/Manhattans/AtCol0_top1SNP_GeneCounts.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=Gene))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Number of Genes")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()
#----------------------------------------------------------------
setwd("~/Projects/BcAt_RNAGWAS")
mydat_c9 <- mydat[mydat$chr==9,]
jpeg("plots/Manhattans/5xRand/AtCol0_chr9_top1SNP.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_c9, aes(x=ps, y=(-log10(p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Distance (Mb)", breaks=c(0,1e+06,2e+06,3e+06),labels=c(0,1,2,3))+
    expand_limits(y=0)
)
dev.off()