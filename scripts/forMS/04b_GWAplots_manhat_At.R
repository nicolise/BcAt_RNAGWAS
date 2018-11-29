#Nicole E Soltis
#07/20/18

#--------------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of Bc reads across ALL GENES 

#from GEMMA_lsm/D_GWAplots

rm(list=ls())
#start with a small file
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At/")
mydat01 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
mydat <- mydat01

#Make plotting variables for snp
mydat_plot <- mydat[order(mydat$chr, mydat$ps),]
mydat_plot$Index.s = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#-------------------------------------------------------------------------------------
#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
myColors <- rep(c("darkgreen","palegreen3"), 9)
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/AtCol0_top1SNP_bySNP.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.s, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list(y="-log10(p)", title=NULL))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    #geom_hline(yintercept=-log10(x), linetype="dashed")+
    expand_limits(y=0)
)
dev.off()

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=FALSE), "Index.s"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=TRUE), "Index.s"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
#-------------------------------------------------------------------------------------
#plot by transcript center!
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(mydat_plot$chr.t)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)
setwd("~/Projects/BcAt_RNAGWAS")
#jpeg("plots/Manhattans/BcLSM_top1SNP_tsdist_byTrans.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.t, y=(ts_dist)))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr.t),alpha=0.001))+
        labs(list(y="Distance (Mb) Transcript Center to Top SNP Hit", title=NULL))+
    scale_y_continuous(breaks= c(5e+05, 1e+06, 1.5e+06, 2e+06, 2.5e+06, 3e+06, 3.5e+06, 4e+06, 4.5e+06), labels=c("0.5", "1", "1.5", "2", "2.5","3","3.5","4","4.5"))+
        theme(legend.position="none")+
        scale_x_continuous(name="Chromosome", breaks = c(2055474, 5781596, 9056392, 11892395, 14584073, 17378738, 20048103, 22665684, 25239020, 27690934, 30056515, 32405141, 34709229, 36904810, 39001984, 40968320, 42053636, 42298445), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))+
        expand_limits(y=0)
  )
dev.off()

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.t, mydat_plot$tmid),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.t, fromLast=FALSE), "Index.t"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.t, fromLast=TRUE), "Index.t"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
#-----------------------------------------------------------------------------------

#manhattan plot: plot SNP location vs. effect estimate for each experiment
#cis vs. trans manhattan plot: plot SNP location vs. distance to gene center for each experiment

#------------------------------------------------------------------------------
#plot Chr 1 only
mydat_plot_c1 <- mydat_plot[mydat_plot$chr.snp==1,]
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_bySNP_Chr1.jpg", width=8, height=5, units='in', res=600)
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

#------------------------------------------------------------------------------------------
#manhattish plot: for each SNP location, # Bc transcripts that have that SNP as the top hit
#new summary df: for each unique value of Index.s, add a variable counting the number of occurrences (= multiple transcripts)
mydat_summ <- mydat_plot[,c("chr","ps","p_score","Gene","Index.s")]
mydat_summ_ngene <- aggregate(Gene ~ Index.s, data = mydat_summ, FUN = function(x){NROW(x)})
#now add SNP data back on, matching by Index.s
mydat_summ_ngene <- merge(mydat_summ_ngene, mydat_summ[,c("chr","ps","Index.s")], by="Index.s")
#remove duplicate rows
mydat_summ_ngene <- unique(mydat_summ_ngene)

library(ggplot2)
#create a custom color scale
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_summ_ngene$chr.snp)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_numgenesPsnp_Bc.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_summ_ngene, aes(x=Index.s, y=(GeneNoTranscript)))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0)+h
    geom_point(aes(color = factor(chr.snp),alpha=0.001))+
    labs(list(y="Number of Genes", title=NULL))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()


#-----------------------------------------------------------
#redraw this, but highlighting top At/ Bc overlap hotspots
setwd("~/Projects/BcAt_RNAGWAS/data")
myhots <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/BcAt_topBcSNPGenes_numTranscripts_funcannot_readin.csv")
names(myhots)[1] <- "Gene"
myhots <- myhots[,c(2,1,3,4,5,6,23,24)]
mydat_summ_ngene$chr.snp <- paste(mydat_summ_ngene$chr, mydat_summ_ngene$ps, sep=".")
myhots_plot <- merge(mydat_summ_ngene, myhots, by="chr.snp", all=TRUE)
myhots_plot$Cat <- ifelse(myhots_plot$Gene.B != 0, (ifelse(myhots_plot$Gene.A != 0, "Both", "B")), "A")
myhots_plot$chr <- as.factor(myhots_plot$chr.x)

#below here: works as cat plot
triColors <- c("darkgreen", "navyblue", "purple")
names(triColors) <- levels(as.factor(myhots_plot$Cat))
tricolScale <- scale_color_manual(name = "Cat",values = triColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_numgenesPsnp_At_tophots.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(myhots_plot, aes(x=Index.s, y=(Gene.x)))+
    theme_bw()+
    tricolScale+
    geom_point(aes(color = factor(Cat)))+
    #scale_fill_manual(values=colScale)+
    #scale_colour_manual(values=c("navyblue","darkgreen","purple"))+
    labs(list(y="Number of Genes", title=NULL))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))
)
dev.off()