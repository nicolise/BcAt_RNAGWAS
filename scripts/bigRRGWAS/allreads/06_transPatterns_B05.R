#06_transPatterns_B05.R
#Nicole E Soltis
#02/22/18

#uses input files from 05_readin_bigRRouts_B05.R
#-----------------------------------------------------------------
rm(list=ls())
library("tidyr"); library("ape"); library("ggplot2")
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
allSNP1 <- read.csv("ChrAll_top1SNPperGene.csv")
allSNP10 <- read.csv("ChrAll_top10SNPperGene.csv")
allSNP <- allSNP10

allSNP <- allSNP[,-c(1)]

#get gff3 file info -- did this on PC
setwd("~/Projects/BcGenome/data/ensembl/BO5.10")
setwd("~/Documents/GitRepos/BcGenome/data/ensembl/BO5.10")
#can extract the files from .gz using 7-zip
#then read using R!

my.gene.list <- as.data.frame(unique(allSNP$gene))
names(my.gene.list)[1] <- "mygenes"
my.gene.list$mygenes <- as.character(my.gene.list$mygenes)
#some lazy regex
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$mygenes, ".1.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".2.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".3.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".4.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".5.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".6.HEM"), "[", 1)

#could do this individually within chromosomes if it's reaaaal slow
my.gff <- read.gff("extracted/Botrytis_cinerea.ASM83294v1.38.gff3", na.strings = c(".", "?"))
#now loop over to keep things
#1:length(my.gene.list$justgenes)

Sys.time()
my.gff.genes <- my.gff[1,]
my.gff.genes$transcript <- NA
my.gff.genes <- my.gff.genes[-c(1),]
#batches 1:1778, 1780:7735, 7737:8196, 8198:end
#broke at y=1779, 7736, 8197
for (y in c(1:1778, 1780:7735, 7737:8196, 8198:length(my.gene.list$justgenes))){
  my.current.gene <- my.gff[grep(my.gene.list[y,2], my.gff$attributes), ]
  my.current.gene$transcript <- my.gene.list[y,2]
  my.gff.genes <- rbind(my.gff.genes, my.current.gene)
}
Sys.time()
#for a first go, just keep the "gene" rows
my.transcript.locs <- my.gff.genes[my.gff.genes$type=="gene", ]
Sys.time()
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout")
write.csv(my.gff.genes, "ChrAll_fullgenesList.csv")
#get average gene length
my.transcript.locs$length <- abs(my.transcript.locs$end - my.transcript.locs$start) 
mean(my.transcript.locs$length) #2489 ... call it 2500
max(my.transcript.locs$length) #longest gene is 21,489 ... big!
#---------------------------------------------------------------------------------
#now: plot for each transcript
#my.full.data includes effect sizes for each transcript across whole genome
#first plot: manhattan of each transcript individually, with a bar for gene location
allSNP <- separate(allSNP, outpt.HEM, into = c("Chrom","Pos") ) 
allSNP$Chr.Pos <- paste(allSNP$Chrom, allSNP$Pos, sep=".")
#get gene midpoints
my.transcript.locs$MidGene <- (my.transcript.locs$end + my.transcript.locs$start)/2
my.transcript <- my.transcript.locs[,c("transcript","MidGene", "seqid")]
names(my.transcript)[3] <- "Gene.Chrom"

#sort dataframe rows in order of Chrom, then Pos
my.plotdat <- allSNP
my.plotdat$Chrom <- as.numeric(my.plotdat$Chrom)
#only keep chromosome 1 SNPs
my.plotdat$Pos <- as.numeric(my.plotdat$Pos)
my.plotdat<- my.plotdat[with(my.plotdat, order(Chrom, Pos)), ]

#now for scatter plot, gene vs. SNPs
names(my.plotdat)[1] <- "SNP.Chrom"
names(my.plotdat)[2] <- "SNP.Pos"
names(my.plotdat)[5] <- "mygenes"
names(my.plotdat)[6] <- "SNP.C.P"

my.plotdat <- merge(my.plotdat, my.gene.list, by = "mygenes")
names(my.plotdat)[7] <- "transcript"

#this can include multiple genes per SNP:
my.plotdat <- merge(my.plotdat, my.transcript, by="transcript")
my.plotdat$SNP.Pos <- as.numeric(my.plotdat$SNP.Pos)
my.plotdat$Gene.Pos <- as.numeric(my.plotdat$MidGene)

#if then to add a column for cis color coding
my.plotdat$Cis <- ifelse(my.plotdat$SNP.Pos<(my.plotdat$Gene.Pos-50000), "trans", ifelse(my.plotdat$SNP.Pos>(my.plotdat$Gene.Pos+50000), "trans", "cis"))

my.plotdat$Cis <- ifelse(my.plotdat$SNP.Chrom == my.plotdat$Gene.Chrom & my.plotdat$SNP.Pos>(my.plotdat$Gene.Pos-50000) & my.plotdat$SNP.Pos<(my.plotdat$Gene.Pos+50000), "cis", "trans")

#----------------------------------------------------------------------------
#now, need to add indexing so that SNPs are sequential along the chromosomes
my.plotdat<- my.plotdat[with(my.plotdat, order(SNP.Chrom, SNP.Pos)), ]
#Add indexing
my.plotdat$SNP.Index = NA
lastbase = 0
for (i in unique(my.plotdat$SNP.Chrom)) {
  print(i)
  if (i==1) {
    my.plotdat[my.plotdat$SNP.Chrom==i, ]$SNP.Index=my.plotdat[my.plotdat$SNP.Chrom==i, ]$SNP.Pos
  }	else {
    lastbase=lastbase+max(subset(my.plotdat,my.plotdat$SNP.Chrom==i-1)$SNP.Pos, 1)
    my.plotdat[my.plotdat$SNP.Chrom==i, ]$SNP.Index=my.plotdat[my.plotdat$SNP.Chrom==i, ]$SNP.Pos+lastbase
  }
}
#sanity check: is each chromosome denoted by a continuously increasing segment of indices?
plot(my.plotdat$SNP.Pos, my.plotdat$SNP.Index)

my.plotdat$Gene.Chrom <- as.numeric(my.plotdat$Gene.Chrom)
my.plotdat<- my.plotdat[with(my.plotdat, order(Gene.Chrom, Gene.Pos)), ]
#Add indexing
my.plotdat$Gene.Index = NA
lastbase = 0
for (i in unique(my.plotdat$Gene.Chrom)) {
  print(i)
  if (i==1) {
    my.plotdat[my.plotdat$Gene.Chrom==i, ]$Gene.Index=my.plotdat[my.plotdat$Gene.Chrom==i, ]$Gene.Pos
  }	else {
    lastbase=lastbase+max(subset(my.plotdat,my.plotdat$Gene.Chrom==i-1)$Gene.Pos, 1)
    my.plotdat[my.plotdat$Gene.Chrom==i, ]$Gene.Index=my.plotdat[my.plotdat$Gene.Chrom==i, ]$Gene.Pos+lastbase
  }
}
#sanity check: is each chromosome denoted by a continuously increasing segment of indices?
plot(my.plotdat$Gene.Pos, my.plotdat$Gene.Index)

#then, manhattan plots: in this case, SNP effect size = y, SNP position = x, color = cis or trans
# for finding the chromosome starts and finish for Chromosome Label Location:
my.plotdat<- my.plotdat[with(my.plotdat, order(SNP.Chrom, SNP.Pos)), ]
my.chroms <- as.data.frame(my.plotdat[!duplicated(my.plotdat$SNP.Chrom, fromLast=FALSE), "SNP.Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- my.plotdat[!duplicated(my.plotdat$SNP.Chrom, fromLast=TRUE), "SNP.Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
#is top SNP on same Chr as gene ?
my.plotdat$cisChrom <- ifelse(my.plotdat$SNP.Chrom == my.plotdat$Gene.Chrom, "same_Chrom", "uniq_Chrom")
#---------------------------------------------------------------------
#now a plot!
#whole genome, every top SNP plotted
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/BO5/WGTrans/MAF20_NA10.top10SNPpGene.ManhattanPlot_cistrans.jpg", width=8, height=5, units='in', res=600)
    ggplot(my.plotdat, aes(x=SNP.Index, y=Estimate))+
         theme_bw()+
         geom_point(aes(color = factor(Cis)), alpha=1/4)+
         labs(list(y=expression(paste("Estimated Effect Size"))))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 6, title=""))+
         #theme(legend.position="none")+
         theme(panel.border = element_blank(), 
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(), 
               axis.line = element_line(colour = "black"))+
         #scale for top 1 SNP
         #scale_x_continuous(name="SNP Chromosome", breaks = c(2043601, 5750286,  9025856, 11865743, 14509489, 17251748, 19895731, 22482767, 25004365, 27459837, 29834912, 32179209, 34464786, 36641858, 38721305, 40710035, 41747570, 41907653), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
        #scale for top 10 SNP
      scale_x_continuous(name="SNP Chromosome", breaks = c(2043507, 5755145, 9036754, 11880741, 14580813, 17404665, 20085410, 22712722, 25269812, 27728242, 30101417, 32447328, 34748576, 36941434, 39023188, 41016119, 42100685, 42303332), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
 dev.off()

jpeg("plots/BO5/WGTrans/Scatter_top10SNPs_WG_cistrans.jpg", width=8, height=5, units='in', res=600)
ggplot(my.plotdat, aes(x=SNP.Index, y=Gene.Index))+
  theme_bw()+
  geom_point(aes(colour=Cis), alpha=1/8)+
  labs(list(y="Center of Gene along B05.10 Genome", x="Position of Top SNPs onB05.10 Genome" ))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

my.plotdat$Dist <- abs(my.plotdat$SNP.Index - my.plotdat$Gene.Index)
jpeg("Top1SNPs_cis_hist.jpg", width=7.5, height=5, units='in', res=600)
hist(my.plotdat$Dist,50, col="slateblue3",xlab="Distance")
dev.off()