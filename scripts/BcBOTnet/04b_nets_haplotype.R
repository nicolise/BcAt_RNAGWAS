#Nicole E Soltis
#08/03/18

#-------------------------------------------------------------------------------------
rm(list=ls())
#Hierarchical cluster analysis
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")

myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"

#use pvclust for clustering with significance/ p values
library(pvclust)

#for BOT:
##originally used *_crop, but fixed SNP pairs for PED in *_recrop
binPEDc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_recrop.ped")
#write.csv(binPEDc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_recrop.csv")
myInfoc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_recrop.info")
row.names(binPEDc12_bot) <- binPEDc12_bot$V2
#clusters_bot <- hclust(dist(binPEDc12_bot[, 7:478]), method="average")

#100 bootstraps takes 6 minutes -- 1000 bootstraps = 1 hour?
#pvclust starts here
botSNP_pvclust <- as.data.frame(t(binPEDc12_bot[,7:478]))
#try adding 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
botSNP_pvclust[nrow(botSNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
botSNP_myvar <- botSNP_pvclust
myrow <- nrow(botSNP_myvar)
for (i in 1:96){
  botSNP_myvar[myrow+1,i] <- ifelse(botSNP_myvar[myrow,i]==2, 1, 2)
}

# #or, option to remove isolates with no variation across SNPs (whole block same haplotype)
# ind <- apply(mysubset, 1, var) == 0
# mysubset <- mysubset[!ind,]
# botSNP_novar <- apply(botSNP_pvclust, 2, var) == 0
# botSNP_myvar <- botSNP_pvclust[,!botSNP_novar]

#test with 10 bootstraps
pvclust_bot <- pvclust(botSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=10)
#run with bootstrapping!
mytime <- Sys.time()
pvclust_bot <- pvclust(botSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#plot clustering
plot(pvclust_bot)
pvrect(pvclust_bot, alpha=0.95)
plot(clusters_bot)

#cluster membership of each isolate
#based on pvclust: 2 or 5 groups
clusterCut2_bot <- cutree(clusters_bot, 2)
binPEDc12_bot$Cluster2 <- clusterCut2_bot
clusterCut5_bot <- cutree(clusters_bot, 5)
binPEDc12_bot$Cluster5 <- clusterCut5_bot

#linear model split by cluster membership
#break them up into their major groups and run a linear model and plot the mean phenotype for the major clusters. Phenotype is pathway expression.

#get pathway-level expression for all bot genes per isolate
#simple linear model: pathExp ~ Isolate + Plant + Isolate/Cluster 
setwd("~/Projects/BcAt_RNAGWAS/")
PVisolates <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
PVbot <- PVisolates[,c(1,2)]
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
BotNets <- MyNets[MyNets$Cluster=="BOT",]
BotPhenos <- PhenosNet[,names(PhenosNet) %in% BotNets$Gene]
BotPhenos <- cbind(PhenosNet$Isolate, BotPhenos)
names(BotPhenos)[1] <- "Isolate"
BotPhenos$mean.Pheno <- rowMeans(BotPhenos[,c(2:8)])
##select rows with cluster membership info
BotClusts <- myclusters[,c(1,2)]
BotPhenosPV <- merge(BotPhenos,BotClusts, by = "Isolate")

#linear model: from pvclust, we will try Cluster2 and Cluster5
#from HClust, Cluster3
mybotmod <- aov(mean.Pheno ~ as.factor(Cluster5), data=BotPhenos)
mybotmod <- aov(mean.Pheno ~ as.factor(bot_pvclust), data=BotPhenosPV)
summary(mybotmod)
TukeyHSD(mybotmod) #3 significantly different from 1 and 2, yeeeey!
#no significance at 2-cluster level
#5-cluster: 5 vs. 1, 2, 3, 4
#pvclust: cluster 2 vs. 1

#plot mean phenotype by cluster
library(ggplot2)
##cluster2 is not informative, cluster3 may be (3 differs from 1&2?)
#p <- ggplot(BotPhenosHC, aes(factor(Cluster3), mean.Pheno))
p <- ggplot(BotPhenosPV, aes(factor(bot_pvclust), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Botrydial_PVcluster.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2) + 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#plot cluster phylogeny
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/Botrydial_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_bot, cex=0.5)
dev.off()

#with PV method
##check output filename
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/Botrydial_IsolatesPVCluster_fix.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_bot)
pvrect(pvclust_bot, alpha=0.95)
dev.off()

#---------------------------------------------------------------------
#for BOA
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
binPEDc1_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa.ped")
write.csv(binPEDc1_boa, "02b_Haploview/binMAF20NA10_chr1_boa.csv")
myInfoc1_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa.info")
row.names(binPEDc1_boa) <- binPEDc1_boa$V2
#v6 is still a dummy phenotype -- not a real SNP
clusters_boa <- hclust(dist(binPEDc1_boa[, 7:length(binPEDc1_boa)]), method="average")

#use pvclust for clustering with significance/ p values
library(pvclust)
#pvclust starts here
boaSNP_pvclust <- as.data.frame(t(binPEDc1_boa[,7:251]))
##one option to remove isolates with no SNP variation
# boaSNP_novar <- apply(boaSNP_pvclust, 2, var) == 0
# boaSNP_myvar <- boaSNP_pvclust[,!boaSNP_novar]
#instead, try adding 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
boaSNP_pvclust[nrow(boaSNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
boaSNP_myvar <- boaSNP_pvclust
myrow <- nrow(boaSNP_myvar)
for (i in 1:96){
  boaSNP_myvar[myrow+1,i] <- ifelse(boaSNP_myvar[myrow,i]==2, 1, 2)
}
#run with 1000 bootstraps
mytime <- Sys.time()
pvclust_boa <- pvclust(boaSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#plot clustering
plot(pvclust_boa)
pvrect(pvclust_boa, alpha=0.95)
plot(clusters_boa)

#cluster membership of each isolate
#could try 2, 3, 4 
plot(clusters_boa)
clusterCut2_boa <- cutree(clusters_boa, 2)
binPEDc1_boa$Cluster2 <- clusterCut2_boa


BoaNets <- MyNets[MyNets$Cluster=="BOA",]
BoaPhenos <- PhenosNet[,names(PhenosNet) %in% BoaNets$Gene]
BoaPhenos <- cbind(PhenosNet$Isolate, BoaPhenos)
names(BoaPhenos)[1] <- "Isolate"
BoaPhenos$mean.Pheno <- rowMeans(BoaPhenos[,c(2:14)])
#BoaClusts <- binPEDc1_boa[,c(252,253,254)]
#BoaClusts[,4] <- rownames(BoaClusts)
#names(BoaClusts)[4] <- "Isolate"
#BoaPhenosHC <- merge(BoaPhenos, BoaClusts, by = "Isolate")
PVboa <- PVisolates[,c(1,3)]
BoaPhenosPV <- merge(BoaPhenos,PVboa, by = "Isolate")

#linear model: from plotting, we will try cluster3
myboamod <- aov(mean.Pheno ~ as.factor(Cluster3), data=BoaPhenos)
summary(myboamod)
TukeyHSD(myboamod) #sig pairings: 2 vs. 1, 3

#plot mean phenotype by cluster
library(ggplot2)
##expect cluster3 to be informative
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/BotcynicAcid_PVclust.jpg", sep=""), width=8, height=6, units='in', res=600)
#p <- ggplot(BoaPhenosHC, aes(factor(Cluster3), mean.Pheno))
p <- ggplot(BoaPhenosPV, aes(factor(boa_pvclust), mean.Pheno))
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#plot cluster phylogeny
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/BotcynicAcid_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_boa, cex=0.5)
dev.off()

#with PV method 
##check output name
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVCluster_fix.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa)
pvrect(pvclust_boa, alpha=0.95)
dev.off()
#--------------------------------------------------------------------
#for Net5
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
binPEDc1_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5.ped")
write.csv(binPEDc1_net5, "02b_Haploview/binMAF20NA10_chr1_net5.csv")
myInfoc1_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5.info")
row.names(binPEDc1_net5) <- binPEDc1_net5$V2
#v6 is still a dummy phenotype -- not a real SNP
clusters_net5 <- hclust(dist(binPEDc1_net5[, 7:length(binPEDc1_net5)]), method="average")
plot(clusters_net5)
clusterCut_net5 <- cutree(clusters_net5, 6)
#cluster membership of each isolate
table(clusterCut_net5, binPEDc1_net5$V2)


#use pvclust for clustering with significance/ p values
library(pvclust)
net5SNP_pvclust <- as.data.frame(t(binPEDc1_net5[,7:length(binPEDc1_net5)]))
# #one option: remove isolates with no variation across SNPs (whole block same haplotype)
# net5SNP_novar <- apply(net5SNP_pvclust, 2, var) == 0
# net5SNP_myvar <- net5SNP_pvclust[,!net5SNP_novar]
#try adding 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
net5SNP_pvclust[nrow(net5SNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
net5SNP_myvar <- net5SNP_pvclust
myrow <- nrow(net5SNP_myvar)
for (i in 1:96){
  net5SNP_myvar[myrow+1,i] <- ifelse(net5SNP_myvar[myrow,i]==2, 1, 2)
}
mytime <- Sys.time()
pvclust_net5 <- pvclust(net5SNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#plot clustering
plot(pvclust_net5)
pvrect(pvclust_net5, alpha=0.95)
plot(clusters_net5)

# #cluster membership of each isolate from Hclust
# #could try 2, 3, 4, 5
# plot(clusters_net5)
# clusterCut2_net5 <- cutree(clusters_net5, 2)
# binPEDc1_net5$Cluster2 <- clusterCut2_net5

#assign isolates to clusters
Net5Nets <- MyNets[MyNets$Cluster=="NET5",]
Net5Phenos <- PhenosNet[,names(PhenosNet) %in% Net5Nets$Gene]
Net5Phenos <- cbind(PhenosNet$Isolate, Net5Phenos)
names(Net5Phenos)[1] <- "Isolate"
Net5Phenos$mean.Pheno <- rowMeans(Net5Phenos[,c(2:11)])
Net5Clusts <- binPEDc1_net5[,c(230:233)]
Net5Clusts[,5] <- rownames(Net5Clusts)
names(Net5Clusts)[5] <- "Isolate"
Net5PhenosHC <- merge(Net5Phenos, Net5Clusts, by = "Isolate")
PVnet5 <- PVisolates[,c(1,4)]
Net5PhenosPV <- merge(Net5Phenos,PVnet5, by = "Isolate")

#linear model: from plotting, we will try cluster3
mynet5mod <- aov(mean.Pheno ~ as.factor(Cluster3), data=Net5Phenos)
summary(mynet5mod)
TukeyHSD(mynet5mod) #sig pairings: 3 vs. 1, 2 

mynet5mod <- aov(mean.Pheno ~ as.factor(Cluster5), data=Net5Phenos)
summary(mynet5mod)
TukeyHSD(mynet5mod) #5 vs. 1, 2, 3, 4

#plot mean phenotype by cluster
library(ggplot2)
##expect cluster3 to be informative
#p <- ggplot(Net5PhenosHC, aes(factor(Cluster3), mean.Pheno))
#p <- ggplot(Net5PhenosHC, aes(factor(Cluster5), mean.Pheno))
p <- ggplot(Net5PhenosPV, aes(factor(net5_pvclust), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Net5_PVclust.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#plot cluster phylogeny
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/Net5_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_net5, cex=0.5)
dev.off()

#with PV method 
##check output name
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/Network5_IsolatesPVCluster_fix.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_net5)
pvrect(pvclust_net5, alpha=0.95)
dev.off()

#plots here
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/Net5_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_net5, cex=0.5)
dev.off()

#----------------------------------------------------------------------------
