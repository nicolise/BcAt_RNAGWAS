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
AtSNP$chr_snp <- paste(AtSNP$chr, AtSNP$ps, sep="_")
BcSNP$chr_snp <- paste(BcSNP$chr, BcSNP$ps, sep="_")
BcSNP2 <- BcSNP[,c(3,5,16,15,17)]
BcSNP2$log10p <- -log10(BcSNP2$p_score)
names(BcSNP2) <- c("chr.B","ps.B","Gene.B","p.B","chr_snp","log10p.B")

AtSNP2 <- AtSNP[,c(4,6,17,16,18)]
AtSNP2$log10p <- -log10(AtSNP2$p_score)
names(AtSNP2) <- c("chr.A","ps.A","Gene.A","p.A","chr_snp","log10p.A")

#here, taking only top p-value (and thus top transcript) per SNP
require(dplyr)
attach(BcSNP2)
BcSNP3 <- as.data.frame(BcSNP2 %>% group_by(chr_snp) %>% filter(log10p.B == max(log10p.B)))

#here, taking only top p-value (and thus top transcript) per SNP
attach(AtSNP2)
AtSNP3 <- as.data.frame(AtSNP2 %>% group_by(chr_snp) %>% filter(log10p.A == max(log10p.A)))

#alternatively, count # genes with each SNP as top 1 hit
AtSNP4 <- aggregate(Gene.A ~ chr_snp, data = AtSNP2, FUN = function(x){NROW(x)})
BcSNP4 <- aggregate(Gene.B ~ chr_snp, data = BcSNP2, FUN = function(x){NROW(x)})

#then merge SNP files together
## choose SNP3 or SNP4 for different plots
BcSNP5 <- BcSNP4
AtSNP5 <- AtSNP4
#first, filter out chr_snp sites if in the hipeaks permutation list
setwd("~/Projects/BcAt_RNAGWAS")
hiAtpeaks<- read.csv("data/GEMMA_eachBc_At/06_GEMMAsumm_RAND/PeaksOver20.csv")
AtSNP5 <- AtSNP5[!AtSNP5$chr_snp %in% hiAtpeaks$chr_ps,]
hiBcpeaks <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm_RAND/PeaksOver5.csv")
BcSNP5 <- BcSNP5[!BcSNP5$chr_snp %in% hiBcpeaks$chr_ps,]
AllSNP <- merge(BcSNP5, AtSNP5, by="chr_snp", all=TRUE)

# ##method for p values plot
# #fill in NA with 0
# AllSNP$log10p.A[is.na(AllSNP$log10p.A)] <- 0
# AllSNP$log10p.B[is.na(AllSNP$log10p.B)] <- 0
# #alt method for B
# TopSNPB <- AllSNP[AllSNP$log10p.B >= 6,]
# TopSNPB <- TopSNPB[TopSNPB$log10p.A==0,]
# TopSNPB <- TopSNPB[order(-TopSNPB$log10p.B),]
# TopSNPB <- TopSNPB[1:20,]
# #alt method for A
# TopSNPA <- AllSNP[AllSNP$log10p.A >= 7,]
# TopSNPA <- TopSNPA[TopSNPA$log10p.B==0,]
# TopSNPA <- TopSNPA[order(-TopSNPA$log10p.A),]
# TopSNPA <- TopSNPA[1:20,]
# #overlap
# TopSNPAll <- AllSNP[AllSNP$log10p.A >= 5,]
# TopSNPAll <- TopSNPAll[TopSNPAll$log10p.B >= 5,]
# TopSNPAll$mean.10p <- rowMeans(TopSNPAll[,c("log10p.A", "log10p.B")])
# TopSNPAll <- TopSNPAll[order(-TopSNPAll$mean.10p),]
# TopSNPAll <- TopSNPAll[1:20,]
# TopSNPAll <- TopSNPAll[,1:11]

##method for transcript # plot
AllSNP$Gene.B[is.na(AllSNP$Gene.B)] <- 0
AllSNP$Gene.A[is.na(AllSNP$Gene.A)] <- 0
#method for B
TopSNPB <- AllSNP[AllSNP$Gene.A==0,]
TopSNPB <- TopSNPB[order(-TopSNPB$Gene.B),]
TopSNPB <- TopSNPB[1:20,]
#method for A
TopSNPA <- AllSNP[AllSNP$Gene.B==0,]
TopSNPA <- TopSNPA[order(-TopSNPA$Gene.A),]
TopSNPA <- TopSNPA[1:20,]
#overlap
TopSNPAll <- AllSNP[AllSNP$Gene.A >= 1,]
TopSNPAll <- TopSNPAll[TopSNPAll$Gene.B >= 1,]
TopSNPAll$mean.gene <- rowMeans(TopSNPAll[,c("Gene.B", "Gene.A")])
TopSNPAll <- TopSNPAll[order(-TopSNPAll$mean.gene),]
TopSNPAll <- TopSNPAll[1:30,1:3]

#plot for both methods
TopSNPAll <- rbind(TopSNPA, TopSNPB, TopSNPAll)
TopSNPAll$Cat <- "A"
TopSNPcolor <- TopSNPAll[,c("chr_snp","Cat")]
plotSNP <- merge(AllSNP,TopSNPcolor,by="chr_snp", all=TRUE)
plotSNP$Cat[is.na(plotSNP$Cat)] <- "B"
plotSNP$Cat <- as.factor(plotSNP$Cat)

library(ggplot2)
myColors <- c("navyblue", "royalblue1")
names(myColors) <- levels(plotSNP$Cat)
colScale <- scale_colour_manual(name = "Cat",values = myColors)


# #draw plot with p values
# plotSNP$mean.p <- ave(plotSNP$log10p.A, plotSNP$log10p.B)
# plotSNP <- plotSNP[order(plotSNP$mean.p),]
# library(ggplot2)
# setwd("~/Projects/BcAt_RNAGWAS")
# jpeg("plots/paper/TopSNPs_interspecific.jpg", width=8, height=5, units='in', res=600)
# print(
# ggplot(plotSNP, aes(x=log10p.B, y=log10p.A))+
#   theme_bw()+
#   colScale+
#   #
#   geom_point(aes(color=factor(Cat), alpha=0.001), stroke=0)+
#   theme(legend.position="none")+
#   scale_y_continuous(name=expression(paste(italic("A. thaliana "),"-log" ["10"],"(p)")))+
#   scale_x_continuous(name=expression(paste(italic("B. cinerea "),"-log" ["10"],"(p)")))
# )
# dev.off()

#log plot
plotSNP$logGenB <- log(plotSNP$Gene.B+1)
plotSNP$logGenA <- log(plotSNP$Gene.A+1)


#plot with gene counts
library(ggplot2)
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/paper/TopSNPs_interspecific_transCount_log_rmPermutSNP_thr.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(plotSNP, aes(x=logGenB, y=logGenA))+
    theme_bw()+
    colScale+
    #
    geom_point(aes(color=factor(Cat), alpha=0.001), stroke=0)+
    theme(legend.position="none")+
    scale_y_continuous(name=expression(paste(italic("A. thaliana "), "log Transcript Count")))+
    scale_x_continuous(name=expression(paste(italic("B. cinerea "),"log Transcript Count")))+
    geom_hline(yintercept=log(150), linetype=3)+
    geom_vline(xintercept=log(20), linetype=3)
)
dev.off()

#write out top SNPs list to annotate
setwd("~/Projects/BcAt_RNAGWAS")
#write.csv(TopSNPAll, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_20topSNPs_pval.csv")
write.csv(TopSNPAll, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_20topSNPs_TranscCount_rmPermutSNP.csv")
write.csv(plotSNP, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_allhotspots_TranscCount_rmPermutSNP.csv")
#---------------------------------------------------------------------------------
#zoom in on genes linked to BcBOA7 hotspot
#chr.snp is 1.39263, which is just upstream of the gene (but this also means it's between BcBOA6 and BcBOA7)

#check which genes are associated!
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

AtBOA <- AtSNP[AtSNP$chr.snp=="1.39263",]
BcBOA <- BcSNP[BcSNP$chr.snp=="1.39263",]
AtBOA <- AtBOA[,-c(3)]
myBOA <- rbind(BcBOA,AtBOA)
write.csv(myBOA, "GEMMA_eachAt_Bc/07_TopSNPs/BcAtgenes_BcBOA7list.csv")

