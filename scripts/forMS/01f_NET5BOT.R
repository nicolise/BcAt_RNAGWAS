#Nicole E Soltis
#08/02/18

#--------------------------------------------------------------------------------
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_bot_named <- read.csv("02b_Haploview/binMAF20NA10_chr12_bot_recrop_named.csv", na.strings=c("","NA"))
mySNP_net5_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_net5_fix_named.csv", na.strings=c("","NA"))
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"
#-----------------------------------------------------------------------------------------------
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
#exploratory analysis: is there a SNP that tags the BOA deletion?
#try this: using SNP state for all isolates in 0 - 20000 bp (first 20kb) of Chr1
  #include BOA deletion as a y variable: s/d (snp vs. deletion) for all isolates
  #model BOA deletion as additive effects of all SNPs
  #does one SNP tag it?