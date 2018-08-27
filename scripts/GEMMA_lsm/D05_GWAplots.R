#Nicole E Soltis
#07/20/18

#--------------------------------------------------------------------------------------
#manhattany plot of SNP location for top 1 SNPs of Bc reads across ALL GENES 
#include Bc lsm vs. Bc permut data!

#see D04_PhenoToGene for how GEMMA phenotypes were renamed to original gene

rm(list=ls())
## choose a matching set of permut SNP samples and Bc-lsm SNP samples
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
#start with a small file
mydat <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")
nameddat <- read.csv("05_GEMMAsumm/GeneNames/AllBcgenes_top10SNP_MAF20NA10_GEMMA_kmat1_Indexed.csv")
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/B05_GEMMA_Bc/Bc_permut/")
nameddat<- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")

#annotate gene center location to each gene
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

#correlation plot: plot SNP location vs. gene center for each experiment
#then, can overlay 2 experiments
#clean up variables for matching
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
#plot by SNP location!
library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(mydat_plot$chr.snp)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)
setwd("~/Projects/BcAt_RNAGWAS")
#jpeg("plots/Manhattans/BcLSM_top1SNP_tsdist_bySNP.jpg", width=8, height=5, units='in', res=600)
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
#-------------------------------------------------------------------------------------
#manhattan plot: plot SNP location vs. effect estimate for each experiment
#cis vs. trans manhattan plot: plot SNP location vs. distance to gene center for each experiment
