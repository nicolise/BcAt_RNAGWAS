#Nicole E Soltis
#08/16/18

##adapted from GEMMA_GWAS/D05_GWAplots_tsDist.R
rm(list=ls())

#annotate gene center location to each gene
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]
my.gtf$midgene <- (my.gtf$V4 + my.gtf$V5)/2

#start: top 10 SNP
## choose a matching set of permut SNP samples and Bc-lsm SNP samples
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
#done: top1SNP
##read in unique files here
mydat10 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top10SNPsample.csv")
mydat1 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top1SNPsample.csv")
mydat100 <- read.csv("05_GEMMAsumm/GeneNames/GEMMA_top100SNPsample.csv")

##careful of which SNP set in use
mydat <- mydat10

#correlation plot: plot SNP location vs. gene center for each experiment
#then, can overlay 2 experiments
#clean up variables for matching
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

#keep only 1 record per gene (ALL exons/ CDS) before merging to SNP effect sizes. 
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
mydat_plot_hist <- mydat_plot[mydat_plot$ts_dist > 0,]

#try it scaled by the length of each chromosome
#chr.s and chr.t must be the same for this dataset. And only goes up to 16 because no SNPs on 17, 18
mydat_plot_hist$Chr.end <- 1
mydat_plot_hist$Chr.start <- 1
for (i in 1:16){
  mydat_plot_hist[mydat_plot_hist$chr.t==i,]$Chr.end <- max(mydat_plot_hist[mydat_plot_hist$chr.t==i,]$tmid)
  mydat_plot_hist[mydat_plot_hist$chr.t==i,]$Chr.start <- min(mydat_plot_hist[mydat_plot_hist$chr.t==i,]$tmid)
}

mydat_plot_hist$intraC.tsdist <- mydat_plot_hist$ts_dist/ (mydat_plot_hist$Chr.end - mydat_plot_hist$Chr.start)

hist(mydat_plot_hist$intraC.tsdist)
hist(mydat_plot_hist$ts_dist)

#for half-page: width = 6.5
#for quarter-page: width = 3.25
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/paper/BotcynicAcid_top10SNP_geneDistHist_sm.jpg", sep=""), width=3.25, height=4, units='in', res=600)
ggplot(data=mydat_plot_hist, aes(mydat_plot_hist$ts_dist)) + 
  geom_histogram(fill="slateblue1", col="black", alpha=0.4, breaks=seq(0, 4e+06, by = 100000), aes(y =..density..))+
  labs(x="Distance (Mb)", y="Frequency")+
  geom_density()+
  scale_x_continuous(breaks = c(0, 5e+5, 1e+6, 1.5e+6, 2e+6, 2.5e+6, 3e+6, 3.5e+6, 4e+6),
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4))+
  theme_bw()
dev.off()


setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/paper/BotcynicAcid_top10SNP_geneDistHist_InChr_sm.jpg", sep=""), width=3.25, height=4, units='in', res=600)
ggplot(data=mydat_plot_hist, aes(mydat_plot_hist$intraC.tsdist)) + 
  geom_histogram(fill="slateblue1", col="black", alpha=0.4, breaks=seq(0, 1, by = 0.05), aes(y=..density..))+
  labs(x="Distance (Fraction of Chromosome)", y="Frequency")+
  geom_density()
dev.off()