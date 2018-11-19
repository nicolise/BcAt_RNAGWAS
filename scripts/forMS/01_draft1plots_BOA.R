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