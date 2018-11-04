#Nicole E Soltis
#11/02/18
#calculate hotspot overlap between Bc and At
#---------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/")
#read in top 1 SNP lists from each
#could also use the SNPannot?
BcSNP  <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
head(BcSNP)
AtSNP <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
head(AtSNP)
#add variable for chr.snp, then merge files on this
AtSNP$chr.snp <- paste(AtSNP$chr, AtSNP$ps, sep=".")
BcSNP$chr.snp <- paste(BcSNP$chr, BcSNP$ps, sep=".")
BcSNP2 <- BcSNP[,c(3,5,16,15,17)]
BcSNP2$log10p <- -log10(BcSNP2$p_score)
names(BcSNP2) <- c("chr.B","ps.B","Gene.B","p.B","chr.snp","log10p.B")
#only keep the highest p gene per SNP
require(dplyr)
attach(BcSNP2)
BcSNP3 <- as.data.frame(BcSNP2 %>% group_by(chr.snp) %>% filter(log10p.B == max(log10p.B)))

AtSNP2 <- AtSNP[,c(4,6,17,16,18)]
AtSNP2$log10p <- -log10(AtSNP2$p_score)
names(AtSNP2) <- c("chr.A","ps.A","Gene.A","p.A","chr.snp","log10p.A")
attach(AtSNP2)
AtSNP3 <- as.data.frame(AtSNP2 %>% group_by(chr.snp) %>% filter(log10p.A == max(log10p.A)))

AllSNP <- merge(BcSNP3, AtSNP3, by="chr.snp", all=TRUE)
#could use all=TRUE ??

#fill in NA with 0
AllSNP$log10p.A[is.na(AllSNP$log10p.A)] <- 0
AllSNP$log10p.B[is.na(AllSNP$log10p.B)] <- 0

#select top 20 SNPs for each category
# TopSNPB <- AllSNP[AllSNP$log10p.A==0,]
# TopSNPB <- TopSNPB[order(-TopSNPB$log10p.B),]
# TopSNPB <- TopSNPB[1:20,]

#alt method for B
TopSNPB <- AllSNP[AllSNP$log10p.B >= 6,]
TopSNPB <- TopSNPB[TopSNPB$log10p.A==0,]
TopSNPB <- TopSNPB[order(-TopSNPB$log10p.B),]
TopSNPB <- TopSNPB[1:20,]

# TopSNPA <- AllSNP[AllSNP$log10p.B==0,]
# TopSNPA <- TopSNPA[order(-TopSNPB$log10p.A),]
# TopSNPA <- TopSNPA[1:20,]

#alt method for A
TopSNPA <- AllSNP[AllSNP$log10p.A >= 7,]
TopSNPA <- TopSNPA[TopSNPA$log10p.B==0,]
TopSNPA <- TopSNPA[order(-TopSNPB$log10p.A),]
TopSNPA <- TopSNPA[1:20,]

#overlap
TopSNPAll <- AllSNP[AllSNP$log10p.A >= 5,]
TopSNPAll <- TopSNPAll[TopSNPAll$log10p.B >= 5,]
TopSNPAll$mean.10p <- ave(TopSNPAll$log10p.A, TopSNPAll$log10p.B)
TopSNPAll <- TopSNPAll[order(-TopSNPAll$mean.10p),]
TopSNPAll <- TopSNPAll[1:20,]

TopSNPAll <- TopSNPAll[,1:11]
TopSNPAll <- rbind(TopSNPA, TopSNPB, TopSNPAll)
TopSNPAll$Cat <- "A"
TopSNPcolor <- TopSNPAll[,c("chr.snp","Cat")]
plotSNP <- merge(AllSNP,TopSNPcolor,by="chr.snp", all=TRUE)
plotSNP$Cat[is.na(plotSNP$Cat)] <- "B"
plotSNP$Cat <- as.factor(plotSNP$Cat)

myColors <- c("navyblue", "royalblue1")
names(myColors) <- levels(plotSNP$Cat)
colScale <- scale_colour_manual(name = "Cat",values = myColors)
plotSNP$mean.p <- ave(plotSNP$log10p.A, plotSNP$log10p.B)
plotSNP <- plotSNP[order(plotSNP$mean.p),]

#draw plot
library(ggplot2)
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/paper/TopSNPs_interspecific.jpg", width=8, height=5, units='in', res=600)
print(
ggplot(plotSNP, aes(x=log10p.B, y=log10p.A))+
  theme_bw()+
  colScale+
  #
  geom_point(aes(color=factor(Cat), alpha=0.001), stroke=0)+
  theme(legend.position="none")+
  scale_y_continuous(name=expression(paste(italic("A. thaliana "),"-log" ["10"],"(p)")))+
  scale_x_continuous(name=expression(paste(italic("B. cinerea "),"-log" ["10"],"(p)")))
)
dev.off()

#write out top SNPs list to annotate
write.csv(TopSNPAll, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_20topSNPs.csv")

#-------------------------------------------------------------------------------------
#annotate gene info to these SNPs
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

my_data <- TopSNPAll
my_data$chr.A[is.na(my_data$chr.A)] <- 0
my_data$chr.B[is.na(my_data$chr.B)] <- 0
my_data$chr <- ifelse(my_data$chr.B == 0, (ifelse(my_data$chr.A == 0, 0, my_data$chr.A)), my_data$chr.B)
my_data$snp <- ifelse(is.na(my_data$ps.A), my_data$ps.B, ifelse(is.na(my_data$ps.B), my_data$ps.A, my_data$ps.B))
#read in single SNP list
my.snps <- my_data
  #calculate gene center
  #calculate distance gene center to SNP
  #add gene with min distance
  #range +-1 kb around each snp: lowrange toprange
  #match snp chromosome.id to gene V1
  #1:18 but have no sig SNPs on chr 17, 18 so actually 1:16
  
  #associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
  for (j in 1:16){
    gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
    my.snp.sub <- my.snps[my.snps$chr==j,]
    #this prevents erroring out if no SNPs within a certain chromosome
    ifelse(nrow(my.snp.sub) == 0, my.snp.sub <- my.snps[1,], my.snp.sub <- my.snp.sub)
    #within these matched sets...
    gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
    for (i in c(1:nrow(my.snp.sub))){
      this.snp <- as.numeric(my.snp.sub[i,"snp"])
      gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
      my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
      this.gene <- gtf.sub[my.closest.gene,]
      #this prevents erroring out if no SNPs within a certain chromosome
      ifelse(nrow(this.gene) == 0, this.gene <- gtf.sub[1,], this.gene <- this.gene)
      this.line <- cbind(my.snp.sub[i,], this.gene)
      this.line$closest.end <- NA
      ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
    }
    ifelse(j == 1, full.snp.genes <- all.genes, full.snp.genes <- rbind(full.snp.genes, all.genes))
    #now only keep genes if nearest end is within +-1kb of SNP (2kb window)
    full.snp.genes$closest.end <- pmin(abs(full.snp.genes$snp - full.snp.genes$V4),abs(full.snp.genes$snp - full.snp.genes$V5)) 
    full.genes.sub <- full.snp.genes[full.snp.genes$closest.end < 1000,]
  }

setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/07_TopSNPs")
write.csv(full.snp.genes, "BcAt_topBcGenes.csv")
#-----------------------------------------------------------------------------------------
