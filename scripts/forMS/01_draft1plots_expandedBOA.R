#Nicole E Soltis
#expanded BOA indel analysis and pvclustering

#-----------------------------------------------------------------------
##color palette:
#figure X1, slateblue1

##Figure X1. BOA violin plot grouped by genic SNPs and deletion info. from script BcBOTnet/08_haplotype_withGenes.R
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"

#now, mean phenotype plots by cluster membership
#assign isolates to clusters
BoaNets <- MyNets[MyNets$Cluster=="BOA",]
BoaPhenos <- PhenosNet[,names(PhenosNet) %in% BoaNets$Gene]
BoaPhenos <- cbind(PhenosNet$Isolate, BoaPhenos)
names(BoaPhenos)[1] <- "Isolate"
BoaPhenos$mean.Pheno <- rowMeans(BoaPhenos[,c(2:length(BoaPhenos))])
#need a new file of cluster membership
BoaClusts <- myclusters[,c("Isolate","boa_pv_gene_indel")] #can also do with boa_pv_gene
BoaPhenosPV <- merge(BoaPhenos,BoaClusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(BoaPhenosPV, aes(factor(boa_pv_gene_indel), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/paper/BotcynicAcid_PVclust_genic_indel.jpg", sep=""), width=6.5, height=4, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="slateblue1", alpha=0.4) + geom_jitter(height = 0, width = 0.1, size=1, alpha=0.5) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2, alpha = 0.8)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#------------------------------------------------------------------
#draw a manhattan plot zooming in on the top SNPs in the BOA region +- 2 kb
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
my.gtf <- my.gtf[,c(1:8,10,12,14)]
my.gtf.c1 <- my.gtf[my.gtf$V1=="Chromosome1",]
my.gtf.boa <- my.gtf[1:158,]

setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
mydat01 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")

mydat_plot <- mydat01
mydat_plot$chr_pos <- paste(mydat_plot$chr, mydat_plot$ps, sep="_")
mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
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
#SNP range out past deletion end: 4029 to 108994 (82614 if just end of cluster)
mydat_plot <- mydat_plot[mydat_plot$chr==1,]
mydat_plot <- mydat_plot[mydat_plot$ps <= 110994,]

library(ggplot2)
setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_bySNP_BOAregion.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    #colScale+
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list(y="-log10(p)", title=NULL))+
    theme(legend.position="none")+
   # scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()
  
#-----------------------------------------------------------------------------
#also try drawing a manhattan plot of the BOA transcript expression within this region?
#BOA = Bcin01g00030.1 to Bcin01g00130.1
#grab these expression profiles from Linux, then plot
  #ideally, range is Chr1, 2029 bp to 110994.
  