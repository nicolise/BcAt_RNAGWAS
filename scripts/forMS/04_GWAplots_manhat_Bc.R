#Nicole E Soltis
#07/20/18

#--------------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of Bc reads across ALL GENES 

#from GEMMA_lsm/D_GWAplots

rm(list=ls())
#start with a small file
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
mydat01 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
#mydat10 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top10SNPsample.csv")
#mydat100 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top100SNPsample.csv")

##check which file I'm reading in!
mydat <- mydat01

#annotate gene center location to each gene
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

#correlation plot: plot SNP location vs. gene center for each experiment
#then, can overlay 2 experiments
#clean up variables for matching
my.gtf$midgene <- (my.gtf$V4 + my.gtf$V5)/2
my.gtf.genes <- my.gtf[,c("V1","V3","V4","V5","V10","V12","V14")]
my.gtf.genes$V12 <- gsub("^gene:","", my.gtf.genes$V12)
my.gtf.genes$V10 <- gsub("^transcript:","", my.gtf.genes$V10)
names(my.gtf.genes) <- c("chr", "exon","start","stop","transcript","Gene","GeneName")
my.gtf.genes$chr <- gsub("^Chromosome","", my.gtf.genes$chr)
#names(mydat.genes)[10] <- "transcript"
#names(mydat.genes)[11] <- "Gene"
library("dplyr")
my.gtf_genesumm <- my.gtf.genes %>%
  group_by(transcript, GeneName, chr) %>%
  summarize(tstart = min(start, na.rm=TRUE),
                         tstop = max(stop, na.rm=TRUE))
my.gtf_genesumm <- as.data.frame(my.gtf_genesumm)
my.gtf_genesumm$tmid <- (my.gtf_genesumm$tstop + my.gtf_genesumm$tstart)/2
names(my.gtf_genesumm)[3] <- "chr.t"

##run this for each SNP input
mydat <- mydat[,-c(1)]
mydat.genes <- mydat[,c(2,4,1,11:ncol(mydat))]
mydat.genes$GeneNoTranscript <- gsub("\\.[0-9]$", '', mydat.genes$Gene)
mydat.genes$transcript <- mydat.genes$Gene
mydat_plot <- merge(mydat.genes, my.gtf_genesumm, by="transcript")

#manhattan plot: location of actual transcript center (phenotype) on x axis, SNP -log10(p) on y axis
#if SNP is on chromosome other than transcript, set distance to 5,000,000: max SNP pos on a chromosome is 4080738 and max gene center is 4104646
#or, trying with 0 for now
names(mydat_plot)[2] <- "chr.snp"

#Make plotting variables for transcript
mydat_plot$chr.t <- as.numeric(mydat_plot$chr.t)
mydat_plot <- mydat_plot[order(mydat_plot$chr.t, mydat_plot$tmid),]
mydat_plot$Index.t = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr.t)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    mydat_plot[mydat_plot$chr.t==i, ]$Index.t=mydat_plot[mydat_plot$chr.t==i, ]$tmid
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr.t==i-1)$tmid, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    mydat_plot[mydat_plot$chr.t==i, ]$Index.t=mydat_plot[mydat_plot$chr.t==i, ]$tmid+lastbase
  }
}

#Make plotting variables for snp
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
mydat_plot$Index.s = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr.snp)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr.snp==i, ]$Index.s=mydat_plot[mydat_plot$chr.snp==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr.s==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr.s==i, ]$Index.s=mydat_plot[mydat_plot$chr.s==i, ]$ps+lastbase
  }
}

#-------------------------------------------------------------------------------------
#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$chr.snp)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_bySNP_thr.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.s, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0)+
    geom_point(aes(color = factor(chr.snp),alpha=0.001))+
    labs(list(y="-log10(p)", title=NULL))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    geom_hline(yintercept=-log10(2.949698e-08), linetype="dashed")+
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
mydat_summ <- mydat_plot[,c("chr.snp","ps","p_score","GeneNoTranscript","Index.s")]
mydat_summ_ngene <- aggregate(GeneNoTranscript ~ Index.s, data = mydat_summ, FUN = function(x){NROW(x)})
#now add SNP data back on, matching by Index.s
mydat_summ_ngene <- merge(mydat_summ_ngene, mydat_summ[,c("chr.snp","ps","Index.s")], by="Index.s")
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
    #, stroke=0)+
    geom_point(aes(color = factor(chr.snp),alpha=0.001))+
    labs(list(y="Number of Genes", title=NULL))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()
