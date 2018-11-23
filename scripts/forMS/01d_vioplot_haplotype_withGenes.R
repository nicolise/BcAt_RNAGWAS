#Nicole E Soltis
#08/02/18

#--------------------------------------------------------------------------------
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_bot_named <- read.csv("02b_Haploview/binMAF20NA10_chr12_bot_recrop_named.csv", na.strings=c("","NA"))
mySNP_boa_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv", na.strings=c("","NA"))
mySNP_net5_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_net5_fix_named.csv", na.strings=c("","NA"))
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"
#-----------------------------------------------------------------------------------
#for Net5
#use pvclust for clustering with significance/ p values
library(pvclust)
#remove SNPs outside of genes
library(dplyr)
net5SNP_pvclust <- mySNP_net5_named %>% 
  select_if(~ !any(is.na(.)))
net5SNP_pvclust <- as.data.frame(t(net5SNP_pvclust[4:nrow(net5SNP_pvclust),]))
names(net5SNP_pvclust) <- unlist(net5SNP_pvclust[1,])
net5SNP_pvclust <- net5SNP_pvclust[-c(1),]
#try add 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
#make pvclust numeric, not a factor
for (i in 1:96){
  net5SNP_pvclust[,i] <- as.numeric(net5SNP_pvclust[,i])
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

#with PV method 
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/Network5_IsolatesPVCluster_onlyGenic.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_net5)
pvrect(pvclust_net5, alpha=0.95)
dev.off()

#now, mean phenotype plots by cluster membership
#assign isolates to clusters
Net5Nets <- MyNets[MyNets$Cluster=="NET5",]
Net5Phenos <- PhenosNet[,names(PhenosNet) %in% Net5Nets$Gene]
Net5Phenos <- cbind(PhenosNet$Isolate, Net5Phenos)
names(Net5Phenos)[1] <- "Isolate"
Net5Phenos$mean.Pheno <- rowMeans(Net5Phenos[,c(2:11)])
#need a new file of cluster membership
Net5Clusts <- myclusters[,c("Isolate","net5_pv_gene")]
Net5PhenosPV <- merge(Net5Phenos,Net5Clusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(Net5PhenosPV, aes(factor(net5_pv_gene), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Net5_PVclust_genic.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()
#----------------------------------------------------------------------
#for Bot
#use pvclust for clustering with significance/ p values
library(pvclust)
#remove SNPs outside of genes
library(dplyr)
botSNP_pvclust <- mySNP_bot_named %>% 
  select_if(~ !any(is.na(.)))
botSNP_pvclust <- as.data.frame(t(botSNP_pvclust[4:nrow(botSNP_pvclust),]))
names(botSNP_pvclust) <- unlist(botSNP_pvclust[1,])
botSNP_pvclust <- botSNP_pvclust[-c(1),]
#try add 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
#make pvclust numeric, not a factor
for (i in 1:96){
  botSNP_pvclust[,i] <- as.numeric(botSNP_pvclust[,i])
}
botSNP_pvclust[nrow(botSNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
botSNP_myvar <- botSNP_pvclust
myrow <- nrow(botSNP_myvar)
for (i in 1:96){
  botSNP_myvar[myrow+1,i] <- ifelse(botSNP_myvar[myrow,i]==2, 1, 2)
}
mytime <- Sys.time()
pvclust_bot <- pvclust(botSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#with PV method 
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/Botrydial_IsolatesPVCluster_onlyGenic.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_bot)
pvrect(pvclust_bot, alpha=0.95)
dev.off()

#now, mean phenotype plots by cluster membership
#assign isolates to clusters
BotNets <- MyNets[MyNets$Cluster=="BOT",]
BotPhenos <- PhenosNet[,names(PhenosNet) %in% BotNets$Gene]
BotPhenos <- cbind(PhenosNet$Isolate, BotPhenos)
names(BotPhenos)[1] <- "Isolate"
BotPhenos$mean.Pheno <- rowMeans(BotPhenos[,c(2:length(BotPhenos))])
#need a new file of cluster membership
BotClusts <- myclusters[,c("Isolate","bot_pv_gene")] 
BotPhenosPV <- merge(BotPhenos,BotClusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(BotPhenosPV, aes(factor(bot_pv_gene), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Bot_PVclust_genic.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#---------------------------------------------------------------------------
#for BoA
#use pvclust for clustering with significance/ p values
library(pvclust)
#remove SNPs outside of genes
library(dplyr)
#clustering includes only genic regions!
boaSNP_pvclust <- mySNP_boa_named %>% 
  select_if(~ !any(is.na(.)))
boaSNP_pvclust <- as.data.frame(t(boaSNP_pvclust[4:nrow(boaSNP_pvclust),]))
names(boaSNP_pvclust) <- unlist(boaSNP_pvclust[1,])
boaSNP_pvclust <- boaSNP_pvclust[-c(1),]
#try add 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
#make pvclust numeric, not a factor
for (i in 1:96){
  boaSNP_pvclust[,i] <- as.numeric(boaSNP_pvclust[,i])
}
boaSNP_pvclust[nrow(boaSNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
boaSNP_myvar <- boaSNP_pvclust
myrow <- nrow(boaSNP_myvar)
for (i in 1:96){
  boaSNP_myvar[myrow+1,i] <- ifelse(boaSNP_myvar[myrow,i]==2, 1, 2)
}
mytime <- Sys.time()
pvclust_boa <- pvclust(boaSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#with PV method 
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVCluster_onlyGenic.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa)
pvrect(pvclust_boa, alpha=0.95)
dev.off()

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
jpeg(paste("plots/Cluster_Means/BotcynicAcid_PVclust_genic_indel.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#-------------------------------------------------------------------------
#now read in without ~incorrect SNP calls. removed 2 SNPs that look like a complete flip from the neighboring regions within gene cluster, 10 SNPs from outside cluster, and removed all nonzero calls from the deletion regions
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_named <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa_rmcalls.csv", na.strings=c("","NA"), header=FALSE)

#for BoA
#use pvclust for clustering with significance/ p values
library(pvclust)
#remove SNPs outside of genes
library(dplyr)
#clustering includes only genic regions! 
#remove column if row "gene" contains NA
boaSNP_pvclust <- mySNP_boa_named[,!is.na(mySNP_boa_named[1,])]
boaSNP_pvclust <- as.data.frame(t(boaSNP_pvclust[4:nrow(boaSNP_pvclust),]))
names(boaSNP_pvclust) <- unlist(boaSNP_pvclust[1,])
boaSNP_pvclust <- boaSNP_pvclust[-c(1),]
#try add 2 random SNPs to prevent isolate dropping
blah <- rbinom(n=96, size=1, prob=0.5)
#this is not a clever way to replace 0 with 2
for (i in 1:96){
  ifelse(blah[i]==0, blah[i] <- 2, blah[i] <- 1)
}
#make pvclust numeric, not a factor
for (i in 1:96){
  boaSNP_pvclust[,i] <- as.numeric(boaSNP_pvclust[,i])
}
boaSNP_pvclust[nrow(boaSNP_pvclust)+1,] <- blah
#add one more row that is the inverse of the row above
boaSNP_myvar <- boaSNP_pvclust
myrow <- nrow(boaSNP_myvar)
for (i in 1:96){
  boaSNP_myvar[myrow+1,i] <- ifelse(boaSNP_myvar[myrow,i]==2, 1, 2)
}
mytime <- Sys.time()
pvclust_boa <- pvclust(boaSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#with PV method 
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVCluster_onlyGenic_rmcalls_clean.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa, print.num=FALSE, print.pv=FALSE)
pvrect(pvclust_boa, alpha=0.95)
dev.off()

#now, mean phenotype plots by cluster membership
#assign isolates to clusters
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"

BoaNets <- MyNets[MyNets$Cluster=="BOA",]
BoaPhenos <- PhenosNet[,names(PhenosNet) %in% BoaNets$Gene]
BoaPhenos <- cbind(PhenosNet$Isolate, BoaPhenos)
names(BoaPhenos)[1] <- "Isolate"
BoaPhenos$mean.Pheno <- rowMeans(BoaPhenos[,c(2:length(BoaPhenos))])
#need a new file of cluster membership
#previously boa_pv_fulldel_onlygene_rmcalls
BoaClusts <- myclusters[,c("Isolate","six_pv_fulldel_onlygene_rmcalls","Iso_labels")] #can also do with boa_pv_gene
BoaPhenosPV <- merge(BoaPhenos,BoaClusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(BoaPhenosPV, aes(factor(six_pv_fulldel_onlygene_rmcalls), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/BotcynicAcid_PVclust_genic_indel_rmcalls_clean.jpg", sep=""), width=6.5, height=4, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="royalblue1", alpha=0.4) + geom_jitter(height = 0, width = 0.1, size=1, alpha=0.5) + 
  geom_text(aes(label=Iso_labels), hjust=-1, vjust=0, size=2, alpha = 0.8, position=position_jitter(width=0.2, height=0))+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#need a new file of ANOVA cluster membership
#includes only the 3 large clusters, excludes low tails of each cluster
AoVclusts <- myclusters[,c("Isolate","boa_anova")] #can also do with boa_pv_gene
BoaPhenosAoV <- merge(BoaPhenos,AoVclusts, by = "Isolate")
fit <- aov(mean.Pheno ~ boa_anova, data=BoaPhenosAoV)
summary(fit) #N.S. p ~ 0.55