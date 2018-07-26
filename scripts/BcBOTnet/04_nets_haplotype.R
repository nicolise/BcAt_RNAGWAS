#Nicole E Soltis
#06/20/18

#----------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreads_bigRR/B05.10/04_NetSubset")
setwd("~/Projects/BcAt_RNAGWAS/data/allreads_bigRR/B05.10/04_NetSubset")
fullfile <- read.csv("BotBoaNet_allGenes_beta.csv")

#java JRE needed for Haploview
#find java installation
#for %i in (java.exe) do @echo.   %~$PATH:i 

#prep haploview input files: split by chromosome

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
# #try just C1 and C12 for BotBoaNet5
# myPED <- read.table("01_PLINK/dpcharMAF20NA10.ped")
# myMAP <- read.table("01_PLINK/dpcharMAF20NA10.map")
# myInfo <- myMAP[,c(2,4)]
# write.table(myInfo, "02b_Haploview/charMAF20NA10.info", col.names=FALSE, row.names=FALSE)
# write.table(myPED, "02b_Haploview/charMAF20NA10.ped", col.names=FALSE, row.names=FALSE)
# #haploview error: >2 alleles at marker 1987 -- use binary coded data!

#try haploview format:
#*.info is 2 tab-separated columns with no headers. C1 is SNP names, C2 is SNP positions
#*.ped is a normal .ped file with no headers. 
#6 columns of labels, then paired SNP states

binPED <- read.table("01_PLINK/dpbinMAF20NA10.ped")
binMAP <- read.table("01_PLINK/dpbinMAF20NA10.map")
#make Chr1 only and Chr12 only mini-files
binMAPc1 <- binMAP[binMAP$V1==1,]
binMAPc12 <- binMAP[binMAP$V1==12,]

#for PED, column 1:6 are phenotypes and column 7:543504 are paired genotypes
binPED[1:10,1:10]
#binMAPc1 is 23350 SNPs which should be the first (23350 * 2) + 6 = 46706 columns of binPED
binPEDc1 <- binPED[,c(1:46706)]
write.table(binPEDc1, "02b_Haploview/binMAF20NA10_chr1.ped", col.names=FALSE, row.names=FALSE)
myInfoc1 <- binMAPc1[,c(2,4)]
write.table(myInfoc1, "02b_Haploview/binMAF20NA10_chr1.info", col.names=FALSE, row.names=FALSE)

binPEDc1 <- read.table("02b_Haploview/binMAF20NA10_chr1.ped")
myInfoc1 <- read.table("02b_Haploview/binMAF20NA10_chr1.info")

#boa limits C1 5,429	61,522
#find which row for each limit in info
#min
myInfoc1.b <- myInfoc1[myInfoc1$V2 < 5429,]
#closest SNP outside region is 5133
#this only works if exact match
which(grepl(5133, myInfoc1$V2)) #row 11
#calculate which column for each in ped
2*11 + 6
#keep column: 1:6, 28:x
#max
myInfoc1.b <- myInfoc1[myInfoc1$V2 > 61522,]
#closest SNP outside region is 61757
#this only works if exact match
which(grepl(61757, myInfoc1$V2)) #row 133
#calculate which column for each in ped
2*133 + 6
#keep column: 1:6, 28:272
binPEDc1_boa <- binPEDc1[,c(1:6, 28:272)]
myInfoc1_boa <- myInfoc1[11:133,]
write.table(binPEDc1_boa, "02b_Haploview/binMAF20NA10_chr1_boa.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_boa, "02b_Haploview/binMAF20NA10_chr1_boa.info", col.names=FALSE, row.names=FALSE)

#net5 limits C1 4,026,579	4,073,855
#find which row for each limit in info
#min
myInfoc1.b <- myInfoc1[myInfoc1$V2 < 4026579,]
mymax <- max(myInfoc1.b$V2)
myminrow <- which(grepl(mymax, myInfoc1$V2))
#calculate which column for each in ped
mymincol <- 2*myminrow + 6
#max
myInfoc1.b <- myInfoc1[myInfoc1$V2 > 4073855,]
mymin <- min(myInfoc1.b$V2)
mymaxrow <- which(grepl(mymin, myInfoc1$V2))
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 6
#keep column: 1:6, 46382:46604
binPEDc1_net5 <- binPEDc1[,c(1:6, 46382:46604)]
myInfoc1_net5 <- myInfoc1[23188:23299,]
write.table(binPEDc1_net5, "02b_Haploview/binMAF20NA10_chr1_net5.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_net5, "02b_Haploview/binMAF20NA10_chr1_net5.info", col.names=FALSE, row.names=FALSE)

#bot limits C12
#for PED, column 1:6 are phenotypes and column 7:543504 are paired genotypes
#get C12 PED limits
binMAP$V5 <- paste(binMAP$V1, binMAP$V4, sep=".")
min(binMAPc12$V4)
which(grepl("12.1156", binMAP$V5))
binMAP[201324,] #yes this one
max(binMAPc12$V4)
which(grepl("12.2351099", binMAP$V5))
binMAP[216362,] #also this one
201324*2 + 6
216362*2 + 6
binPEDc12 <- binPED[,c(1:6, 402654:432730)]
myInfoc12 <- binMAPc12[,c(2,4)]
write.table(binPEDc12, "02b_Haploview/binMAF20NA10_chr12.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc12, "02b_Haploview/binMAF20NA10_chr12.info", col.names=FALSE, row.names=FALSE)

binPEDc12 <- read.table("02b_Haploview/binMAF20NA10_chr12.ped")
myInfoc12 <- read.table("02b_Haploview/binMAF20NA10_chr12.info")

#bot limits C12 2,217,306	2,245,431
#min
myInfoc12.b <- myInfoc12[myInfoc12$V4 < 2245431,]
mymax <- max(myInfoc12.b$V4)
myminrow <- which(grepl(mymax, myInfoc12$V4))
#calculate which column for each in ped
mymincol <- 2*myminrow + 5
#max
myInfoc12.b <- myInfoc12[myInfoc12$V4 > 2217306,]
mymin <- min(myInfoc12.b$V4)
mymaxrow <- which(grepl(mymin, myInfoc12$V4))
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 5
#keep column: 1:6, 28702:28228
binPEDc12_bot <- binPEDc12[,c(1:6, 28228:28702)]
myInfoc12_bot <- myInfoc12[c(14111:14348),]
#remove last 2 snps from each for haploview - appear to be outside high LD region
#rm snp 2243688
binPEDc12_bot <- binPEDc12[,c(1:6, 28228:28699)]
myInfoc12_bot <- myInfoc12[c(14111:14345),]

write.table(binPEDc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_crop.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_crop.info", col.names=FALSE, row.names=FALSE)
#-----------------------------------------------------------------------------------
#gene clusters +- 1 kb -- for BOA LD does not break down by 1kb, so going out to 2kb
#boa limits C1 5,429	61,522
#find which row for each limit in info
#min
#minimum is 4029 which is less than 2kb out
myInfoc1.b <- myInfoc1[myInfoc1$V2 == 4029,]
mymax <- max(myInfoc1.b$V2)
myminrow <- which(grepl(mymax, myInfoc1$V2))
myminrow <- 1
#calculate which column for each in ped
mymincol <- 2*myminrow + 5
#max
myInfoc1.b <- myInfoc1[myInfoc1$V2 > 63522,]
mymin <- min(myInfoc1.b$V2)
mymaxrow <- which(grepl(mymin, myInfoc1$V2))
mymaxrow <- 158
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 5
#keep column: 1:6, 22:312
binPEDc1_boa_2kb <- binPEDc1[,c(1:6, 7:321)]
myInfoc1_boa_2kb <- myInfoc1[1:158,]
write.table(binPEDc1_boa_2kb, "02b_Haploview/binMAF20NA10_chr1_boa_2kb.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_boa_2kb, "02b_Haploview/binMAF20NA10_chr1_boa_2kb.info", col.names=FALSE, row.names=FALSE)

#net5 limits C1 4,026,579	4,073,855
#find which row for each limit in info
#min
myInfoc1.b <- myInfoc1[myInfoc1$V2 < 4025579,]
mymax <- max(myInfoc1.b$V2)
myminrow <- which(grepl(mymax, myInfoc1$V2))
#calculate which column for each in ped
mymincol <- 2*myminrow + 6
#max
myInfoc1.b <- myInfoc1[myInfoc1$V2 > 4074855,]
mymin <- min(myInfoc1.b$V2)
mymaxrow <- which(grepl(mymin, myInfoc1$V2))
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 6
#keep column: 1:6, 46382:46604
binPEDc1_net5_1kb <- binPEDc1[,c(1:6, 46364:46612)]
myInfoc1_net5_1kb <- myInfoc1[23179:23303,]
write.table(binPEDc1_net5_1kb, "02b_Haploview/binMAF20NA10_chr1_net5_1kb.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_net5_1kb, "02b_Haploview/binMAF20NA10_chr1_net5_1kb.info", col.names=FALSE, row.names=FALSE)

##start here
#bot limits C12 2,217,306	2,245,431
#C12 end is 2,351,099
#min
myInfoc12.b <- myInfoc12[myInfoc12$V2 < 2246431,]
mymax <- max(myInfoc12.b$V2)
myminrow <- which(grepl(mymax, myInfoc12$V2))
#calculate which column for each in ped
mymincol <- 2*myminrow + 6
#max
myInfoc12.b <- myInfoc12[myInfoc12$V2 > 2216306,]
mymin <- min(myInfoc12.b$V2)
mymaxrow <- which(grepl(mymin, myInfoc12$V2))
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 6
#keep column: 1:6, 28186:28724
binPEDc12_bot_1kb <- binPEDc12[,c(1:6, 28186:28724)]
myInfoc12_bot_1kb <- myInfoc12[c(14090:14359),]

write.table(binPEDc12_bot_1kb, "02b_Haploview/binMAF20NA10_chr12_bot_1kb.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc12_bot_1kb, "02b_Haploview/binMAF20NA10_chr12_bot_1kb.info", col.names=FALSE, row.names=FALSE)
#---------------------------------------------------------------------
#Hierarchical cluster analysis
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")

#use pvclust for clustering with significance/ p values
library(pvclust)

#for BOT:
binPEDc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_crop.ped")
myInfoc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_crop.info")
row.names(binPEDc12_bot) <- binPEDc12_bot$V2
clusters_bot <- hclust(dist(binPEDc12_bot[, 7:478]), method="average")

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
BotClusts <- binPEDc12_bot[,c(480,481,482)]
BotClusts[,4] <- rownames(BotClusts)
names(BotClusts)[4] <- "Isolate"
BotPhenosHC <- merge(BotPhenos, BotClusts, by = "Isolate")
BotPhenosPV <- merge(BotPhenos,PVbot, by = "Isolate")

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
p <- ggplot(BotPhenos, aes(factor(Cluster3), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Botrydial_3Cluster.jpg", sep=""), width=8, height=6, units='in', res=600)
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
jpeg(paste("plots/HClust/PVclusters/Botrydial_IsolatesPVCluster_dropIso.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_bot)
pvrect(pvclust_bot, alpha=0.95)
dev.off()

#---------------------------------------------------------------------
#for BOA
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
binPEDc1_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa.ped")
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
clusterCut3_boa <- cutree(clusters_boa, 3)
binPEDc1_boa$Cluster3 <- clusterCut3_boa
clusterCut4_boa <- cutree(clusters_boa, 4)
binPEDc1_boa$Cluster4 <- clusterCut4_boa

BoaNets <- MyNets[MyNets$Cluster=="BOA",]
BoaPhenos <- PhenosNet[,names(PhenosNet) %in% BoaNets$Gene]
BoaPhenos <- cbind(PhenosNet$Isolate, BoaPhenos)
names(BoaPhenos)[1] <- "Isolate"
BoaClusts <- binPEDc1_boa[,c(252,253,254)]
BoaClusts[,4] <- rownames(BoaClusts)
names(BoaClusts)[4] <- "Isolate"
BoaPhenos <- merge(BoaPhenos, BoaClusts, by = "Isolate")
BoaPhenos$mean.Pheno <- rowMeans(BoaPhenos[,c(2:14)])

#linear model: from plotting, we will try cluster3
myboamod <- aov(mean.Pheno ~ as.factor(Cluster3), data=BoaPhenos)
summary(myboamod)
TukeyHSD(myboamod) #sig pairings: 2 vs. 1, 3

#plot mean phenotype by cluster
library(ggplot2)
##expect cluster3 to be informative
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/BotcynicAcid_3Cluster.jpg", sep=""), width=8, height=6, units='in', res=600)
p <- ggplot(BoaPhenos, aes(factor(Cluster3), mean.Pheno))
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
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVCluster_AllIsos.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa)
pvrect(pvclust_boa, alpha=0.95)
dev.off()
#--------------------------------------------------------------------
#for Net5
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
binPEDc1_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5.ped")
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
net5SNP_pvclust <- as.data.frame(t(binPEDc1_net5[,7:233]))
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
Net5Clusts <- binPEDc1_net5[,c(230:233)]
Net5Clusts[,5] <- rownames(Net5Clusts)
names(Net5Clusts)[5] <- "Isolate"
Net5Phenos <- merge(Net5Phenos, Net5Clusts, by = "Isolate")
Net5Phenos$mean.Pheno <- rowMeans(Net5Phenos[,c(2:11)])

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
p <- ggplot(Net5Phenos, aes(factor(Cluster3), mean.Pheno))
p <- ggplot(Net5Phenos, aes(factor(Cluster5), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/Net5_3Cluster.jpg", sep=""), width=8, height=6, units='in', res=600)
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
jpeg(paste("plots/HClust/PVclusters/Network5_IsolatesPVCluster_AllIsos.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_net5)
pvrect(pvclust_net5, alpha=0.95)
dev.off()

#plots here
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/Net5_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_net5, cex=0.5)
dev.off()

#----------------------------------------------------------------------------
#haploview command line notes

#chromosome 1 is too much data for the haploview GUI. can try from command line
#get working directory
#DIR 
#cd to the top of the directory tree (C:)
#cd\ 
#cd "Program Files (x86)\HaploView"
#open the gui
#java .jar HaploView.jar"
#open without gui, get help commands
#java .jar "HaploView.jar" -nogui -help 
#trying -png option: 
#fatal error: Exception in thread "main" java.lang.OutOfMemoryError: java heap space
#trying -blockoutput option:
#fatal error: Exception in thread "main" java.lang.OutOfMemoryError: java heap space
#adding option -mem 1000 (default is 512, -mem 2000 gives an error "could not reserve enough space") 
#still runs out of space with -mem 1000, 1500
#trying with only -blockoutput SPI (rather than -blockoutput ALL)

#java -jar Haploview.jar -nogui -memory 1500 -pedfile NESfiles\binMAF20NA10_chr1.ped -info NESfiles\binMAF20NA10_chr1.info -blockoutput SPI
#----------------------------------------------------------------------------
#trying fastPHASE
#from Gautier 2012
#prep fastPHASE input file
#then run fastPHASE
#then run rehh if I want to look for signatures of selection
library("rehh")
??rehh

#using same genotype input from B05.10 BcAtGWAS bigRR and GEMMA
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#use same SNP set from B05.10 bigRR
mySNPs <- read.csv("B05_GEMMA/01_PLINK/OriginalSNPdata.csv")
#try one for Chr1 only

# #file format for fastPHASE 
# 3
# 4
# # id 1
# 1a11
# 0t01
# # id 2
# 1t11
# 0a00
# # id 3
# ?a01
# ?t10

#line 1: no. individuals
#line 2: no. SNP sites
#line 3: individual id
#line 4: all SNPs, individual 1
#line 5: all SNPs, individual 2 (it's diploid)

myFast <- mySNPs[,4:99]
#transpose and format for PED
myFast2 <- as.data.frame(t(myFast))

#no. individuals = 96
#no. SNP sites = 271749 for WG
  
#this is super fast:
myFast3 <- myFast2[rep(1:nrow(myFast2),each=2),] 

write.table(myFast3, "BcBotGWAS/05_haplotypes/WGS_B05.10_fastPhase.txt", col.names=FALSE)


  
#running fastPHASE
#in Documents/
#./fastPHASE -h (for help options)
# ./fastPHASE -oChr1Try1