#Nicole E Soltis
#07/20/18

#--------------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of Bc reads across ALL GENES 
#include Bc lsm vs. Bc permut data!

#see D04_PhenoToGene for how GEMMA phenotypes were renamed to original gene

rm(list=ls())
## choose a matching set of permut SNP samples and Bc-lsm SNP samples
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
#done: top1SNP
##read in unique files here
mydat100 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top100SNPsample.csv")
mydat10 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top10SNPsample.csv")
mydat1 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")
mydat<- mydat1
#annotate gene center location to each gene
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

my.gtf$midgene <- (my.gtf$V4 + my.gtf$V5)/2
mydat.genes <- mydat[,c(3,5,2,10:ncol(mydat))]
mydat.genes$GeneNoTranscript <- gsub("\\.[0-9]$", '', mydat.genes$Gene)
my.gtf.genes <- my.gtf[,c("V1","V3","V4","V5","V10","V12","V14")]
my.gtf.genes$V12 <- gsub("^gene:","", my.gtf.genes$V12)
my.gtf.genes$V10 <- gsub("^transcript:","", my.gtf.genes$V10)
names(my.gtf.genes) <- c("chr", "exon","start","stop","transcript","Gene","GeneName")
my.gtf.genes$chr <- gsub("^Chromosome","", my.gtf.genes$chr)
names(mydat.genes)[10] <- "transcript"
names(mydat.genes)[11] <- "Gene"
library("dplyr")
my.gtf_genesumm <- my.gtf.genes %>%
  group_by(transcript, GeneName, chr) %>%
  summarize(tstart = min(start, na.rm=TRUE),
                         tstop = max(stop, na.rm=TRUE))
my.gtf_genesumm <- as.data.frame(my.gtf_genesumm)
my.gtf_genesumm$tmid <- (my.gtf_genesumm$tstop + my.gtf_genesumm$tstart)/2
names(my.gtf_genesumm)[3] <- "chr.t"


mydat_plot <- merge(mydat.genes, my.gtf_genesumm, by="transcript")

#manhattan plot: location of actual transcript center on x axis, distance from gene to SNP on y axis
#if SNP is on chromosome other than transcript, set distance to 5,000,000: max SNP pos on a chromosome is 4080738 and max gene center is 4104646
#or, trying with 0 for now
names(mydat_plot)[2] <- "chr.snp"
mydat_plot$ts_dist <- ifelse(mydat_plot$chr.snp == mydat_plot$chr.t, abs(mydat_plot$tmid - mydat_plot$ps), 0)

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

#---------------------------------------------------------------------------------------
#correlation plot: plot SNP location vs. gene center for each experiment
#scatterplots / cis diagonal: transcript center to top SNP!
#options(device = "RStudioGD")
#options(device = "windows") #if want to plot to external window
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(mydat_plot$chr.snp)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#make max distance from gene to snp = 1, then increase from there 
##get max ts_dist
#top 10 = 3791215 
#top 1 SNP = 3671739 
mydat_plot$dist.alpha <- ifelse(mydat_plot$ts_dist == 0, 0, max(mydat_plot$ts_dist) - mydat_plot$ts_dist)

mydat_test <- mydat_plot[1:100,]

setwd("~/Projects/BcAt_RNAGWAS/plots")
##check plot name
jpeg("cisDiagonal/Bclsm_top1SNP_cisdiag.jpg", width=6.5, height=4, units='in', res=600)
print(
ggplot(mydat_plot, aes(x=Index.s, y=Index.t))+
  theme_bw()+
  #size = 0.1, stroke = 0 is too fine for 10 SNP and 1 SNP.
  #for 10 SNP, size = 0.5
  #for 1 SNP, size = 1
  geom_point(aes(color = factor(chr.t),alpha=dist.alpha), size=1, stroke=0)+
  #geom_point(aes(color = factor(chr.t),alpha=dist.alpha), size=0.1)+
  guides(fill=FALSE, color=FALSE)+
  labs(x = "SNP Genomic Position (Mb)", y = "Gene Center Genomic Position (Mb)", alpha="Distance SNP to Gene")+
  scale_x_continuous(breaks = c(0, 1e7, 2e7, 3e7, 4e7), labels=c(0, 10, 20, 30, 40))+
  scale_y_continuous(breaks = c(0, 1e7, 2e7, 3e7, 4e7), labels=c(0, 10, 20, 30, 40))+
  ##check this scale labeling from alpha defaults and alpha = ts_dist
  #for 100 SNP, alpha = 0.05 - 1
  #for 10, 1 SNP = 0.05 - 0.5
  scale_alpha_continuous(range=c(0.05,0.5), breaks=c(0,1e6, 2e6, 3e6), labels=c("Off Chromosome", "3 Mb", "2 Mb", "1 Mb"))
)
dev.off()

#----------------------------------------------------------------------------------
#plot by transcript center!
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(mydat_plot$chr.t)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)
setwd("~/Projects/BcAt_RNAGWAS")
##uniquely name plot here
#jpeg("plots/Manhattans/BcLSM_top10SNP_tsdist_byTrans.jpg", width=8, height=5, units='in', res=600)
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
#plot by SNP location!
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(mydat_plot$chr.snp)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)
setwd("~/Projects/BcAt_RNAGWAS")
##uniquely name plot here
#jpeg("plots/Manhattans/BcLSM_top10SNP_tsdist_bySNP.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.s, y=(ts_dist)))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(chr.snp),alpha=0.001))+
    labs(list(y="Distance (Mb) Transcript Center to Top SNP Hit", title=NULL))+
    scale_y_continuous(breaks= c(5e+05, 1e+06, 1.5e+06, 2e+06, 2.5e+06, 3e+06, 3.5e+06, 4e+06, 4.5e+06), labels=c("0.5", "1", "1.5", "2", "2.5","3","3.5","4","4.5"))+
    theme(legend.position="none")+
    scale_x_continuous(name="Chromosome", breaks = c(2051481, 5749572, 8991404, 11798707, 14515550, 17296053, 19947501, 22566871, 25142418, 27612478, 29978265, 32290406, 34575912, 36762846, 38845087, 40005713), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    expand_limits(y=0)
)
dev.off()

#get chromosome midpoints
mydat_plot <- mydat_plot[order(mydat_plot$chr.snp, mydat_plot$ps),]
my.chroms <- as.data.frame(mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=FALSE), "Index.s"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- mydat_plot[!duplicated(mydat_plot$chr.snp, fromLast=TRUE), "Index.s"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
#---------------------------------------------------------------------------
#what about for SNPs at edges of chromosomes?
head(mydat_plot)
#try: scale each SNP to relative distance from end of chromosome
mydat_plot$Chr.end <- 1
mydat_plot$Chr.start <- 1
for (i in 1:16){
  mydat_plot[mydat_plot$chr.snp==i,]$Chr.end <- max(mydat_plot[mydat_plot$chr.snp==i,]$Index.s)
  mydat_plot[mydat_plot$chr.snp==i,]$Chr.start <- min(mydat_plot[mydat_plot$chr.snp==i,]$Index.s)
}

#are genes only annotated for beginning of chromosomes > 1??
mydat_plot$intraC.sdist <- mydat_plot$ps / (mydat_plot$Chr.end - mydat_plot$Chr.start) #on chr 18, gene centers go past SNP list
ggplot(mydat_plot, aes(x=intraC.sdist, y=(ts_dist)))+
  theme_bw()+
  geom_point(aes(color = factor(chr.t),alpha=0.001))
#more pronounced for SNPs: distance from top SNP to gene appears elevated on average for SNPs at edges of chromosomes

#plot: ts_dist_dc as fraction of max possible ts_dist based on SNP location on chromosome!
#scale ts_dist as relative to chr length
mydat_plot$Chr.end.t <- 1
mydat_plot$Chr.start.t <- 1
#using chromosome delineations from SNPs
for (i in 1:16){
  mydat_plot[mydat_plot$chr.t==i,]$Chr.end.t <- max(mydat_plot[mydat_plot$chr.t==i,]$Index.t)
  mydat_plot[mydat_plot$chr.t==i,]$Chr.start.t <- min(mydat_plot[mydat_plot$chr.t==i,]$Index.t)
}
##does this make any sense?
mydat_plot$intraC.tdist <- mydat_plot$tmid / (mydat_plot$Chr.end.t - mydat_plot$Chr.start.t)
mydat_plot$intraC.tsdist <- abs(mydat_plot$intraC.sdist - mydat_plot$intraC.tdist)
mydat_plot_scale <- mydat_plot[mydat_plot$intraC.sdist <= 1,]
mydat_plot_scale <- mydat_plot_scale[mydat_plot_scale$intraC.tdist <= 1,]
#furthest any location on the chromosome could be from the current SNP
mydat_plot_scale$ts_dist_max <- pmax(abs(mydat_plot_scale$intraC.sdist-1),abs(1-mydat_plot_scale$intraC.sdist))
mydat_plot_scale$ts_dist_dc <- mydat_plot_scale$intraC.tsdist/ mydat_plot_scale$ts_dist_max

ggplot(mydat_plot_scale, aes(x=intraC.tdist, y=(ts_dist_dc)))+
  theme_bw()+
  geom_point(aes(color = factor(chr.t),alpha=0.001))

#-------------------------------------------------------------------------------------
#looking into: transcript density across chromosome, SNP density across chromosome
head(my.gtf_genesumm)
max(my.gtf_genesumm[my.gtf_genesumm$chr.t==1,]$tstop)
max(mydat_plot$Index.t) #42372315
head(mydat_plot)
mydat_plot_gs <- mydat_plot %>%
  group_by(transcript, Gene, GeneName, chr.t) %>%
  summarize(tstart = min(tstart, na.rm=TRUE),
            tstop = max(tstop, na.rm=TRUE),
            tmid = mean(tmid, na.rm=TRUE),
            Index.t = mean(Index.t, na.rm=TRUE),
            ts_dist = mean(ts_dist, na.rm=TRUE))
mydat_plot_gs <- as.data.frame(mydat_plot_gs)

#try: scale each gene to relative distance from end of chromosome
mydat_plot_gs$Chr.end <- 1
mydat_plot_gs$Chr.start <- 1
for (i in 1:18){
  mydat_plot_gs[mydat_plot_gs$chr.t==i,]$Chr.end <- max(mydat_plot_gs[mydat_plot_gs$chr.t==i,]$Index.t)
  mydat_plot_gs[mydat_plot_gs$chr.t==i,]$Chr.start <- min(mydat_plot_gs[mydat_plot_gs$chr.t==i,]$Index.t)
}

mydat_plot_gs$intraC.tdist <- mydat_plot_gs$tmid / (mydat_plot_gs$Chr.end - mydat_plot_gs$Chr.start) #on chr 18, gene centers go past SNP list
#remove C17, 18: no SNPs
mydat_plot_gs <- mydat_plot_gs[!mydat_plot_gs$chr.t==17,]
mydat_plot_gs <- mydat_plot_gs[!mydat_plot_gs$chr.t==18,]
ggplot(mydat_plot_gs, aes(x=intraC.tdist, y=(ts_dist)))+
  theme_bw()+
  geom_point(aes(color = factor(chr.t),alpha=0.001))

#distance from top SNP to gene appears *slightly* elevated on average for transcripts at edges of chromosomes
#---------------------------------------------------------------------------
#try as bar plot/ histogram
ggplot(mydat_plot_gs, aes(Index.t)) + geom_histogram(binwidth = 200) #not informative as is

#try: 1kb windows along Index.t
#if tmid for a line is in that 1kb window -> new column = 1
my_win <- as.data.frame(seq(1,42373000, by=1000))
names(my_win)[1] <- "winStart"
my_win$winStop <- seq(1000, 42373000, by=1000)
#for each row in my_win, if there is an mydat_plot_gs$Index.t between winStart and winStop, variable "Index" = 1. Else = 0.
#for (i in 1:nrow(my_win)){
for (i in 1:10){
  ifelse(mydat_plot_gs$Index.t < my_win$winStop && mydat_plot_gs$Index.t > )  
}