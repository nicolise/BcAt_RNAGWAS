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
#get just BOA gene regions for plot part C
#end line 158
rm(list=ls())
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
my.gtf <- my.gtf[,c(1:8,10,12,14)]
my.gtf.c1 <- my.gtf[my.gtf$V1=="Chromosome1",]
my.gtf.boa <- my.gtf[1:158,]

#SNP range: 4029 to 108994

#draw gene regions as bars -- draw exons
library("genemodel")
data("AT5G62640")
genemodel.plot(model=AT5G62640, start=25149433, bpstop=25152541, orientation="reverse", xaxis=T)
mutation.plot(25150214, 25150218, text="P->S", col="black", drop=-.15, haplotypes=c("red", "blue"))
for (i in c(25150808, 25151000, 25251400)){
  mutation.plot(i, i, col="black", drop=-0.15, haplotypes="red")
}

#for each gene, make a 2-column dataframe with col 1 = type and col 2 = coordinates (range)
my.boa.mods <- my.gtf.boa[,c("V12","V3","V4","V5")]
my.boa.mods$coordinates <- paste(my.boa.mods$V4, my.boa.mods$V5, sep="-")
names(my.boa.mods)[2] <- "type"
my.boa.mods <- my.boa.mods[,c(1,2,5)]
#drop empty levels of gene before splitting df
my.boa.mods$V12 <- droplevels(my.boa.mods$V12)
#replace "CDS" with "coding_region"
split.df <- split(my.boa.mods, my.boa.mods$V12)
#make separate, labeled dfs for each gene
for (i in 1:length(split.df)){
  blah <- unique(split.df[[i]]$V12)
  split.df[[i]] <- split.df[[i]][,c(2:3)]
  assign(paste(blah), split.df[[i]])
}

#messy test version in excel
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
#write.csv(`gene:Bcin01g00010`, "02_BOAnet/test_Bcin01g00010.csv")
Bcin01g00010 <- read.csv("02_BOAnet/test_Bcin01g00010.csv")

genemodel.plot(model=`Bcin01g00010`, start=5429, bpstop=6804, orientation="forward",xaxis=T)

# myColors <- 
# my.d.plot$t <- droplevels(my.d.plot$t)
# names(myColors) <- levels(my.d.plot$t)
# colScale <- scale_colour_manual(name = "gene",values = myColors)
# 
# mycol <- ifelse(col, 'white', 'gray53')
# 
# my.rect$Cols <- NULL
# levels(my.rect$Cols) <- c()
#   
# for (j in 1:length(unique(my.rect$V10))){
#   if((j %% 2) ==0){
#     my.rect$Cols <- "navyblue"
#   } else {
#     my.rect$Cols <- "royalblue1"
#   }
# }
#                  rep(c("navyblue", "royalblue1"),19)
#                  
#                  if((num %% 2) == 0) {
#                    print(paste(num,"is Even"))
#                  } else {
#                    print(paste(num,"is Odd"))
#                  }

setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_named <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa_rmcalls.csv", na.strings=c("","NA"))

snplist <- as.numeric(as.character(t(as.vector(mySNP_boa_named[2,2:length(mySNP_boa_named)]))))

my.rect <- my.gtf.boa[,c(4,5,7,3,9,11)]
y1 <- rep(1, 158)
y2 <- rep(2, 158)
x1 <- my.rect[,"V4"]
x2 <- my.rect[,"V5"]
t <- my.rect[,"V10"]
my.d.plot <- data.frame(x1,x2,y1,y2,t)
plot1 <- ggplot()+
  scale_y_continuous(limits=c(1,20))+
  scale_x_continuous(name="Distance (kb)", breaks=c(25000, 50000, 75000, 100000), labels=c(25, 50, 75, 100))+
  theme_bw()+
  guides(fill=FALSE)+
  colScale+
  geom_rect(data=my.d.plot, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.3)+
  geom_point(aes(x=82614, y=5))+
  geom_point(aes(x=snplist), y=2.5, alpha=0.3)
#add a blip at end of deletion: 82614
#deletion length
82614-4029

#deletion from START to between 160 and 170. 170 is third from last gene
setwd("~/Projects/BcAt_RNAGWAS/plots/paper")
jpeg("BOAgenemodels.jpg", width=8, height=4, units='in', res=600)
plot1
dev.off()

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
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVCluster_onlyGenic_rmcalls.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa)
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
BoaClusts <- myclusters[,c("Isolate","boa_pv_fulldel_onlygene_rmcalls")] #can also do with boa_pv_gene
BoaPhenosPV <- merge(BoaPhenos,BoaClusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(BoaPhenosPV, aes(factor(boa_pv_fulldel_onlygene_rmcalls), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/BotcynicAcid_PVclust_genic_indel_rmcalls.jpg", sep=""), width=6.5, height=4, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="slateblue1", alpha=0.4) + geom_jitter(height = 0, width = 0.1, size=1, alpha=0.5) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2, alpha = 0.8)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()

#need a new file of ANOVA cluster membership
#includes only the 3 large clusters (4,5,1), excludes low tails of each cluster
#boa_anova is group membership, only includes isolates for ANOVA
AoVclusts <- myclusters[,c("Isolate","boa_anova")] #can also do with boa_pv_gene
BoaPhenosAoV <- merge(BoaPhenos,AoVclusts, by = "Isolate")
fit <- aov(mean.Pheno ~ boa_anova, data=BoaPhenosAoV)
summary(fit) #N.S. p ~ 0.55

#-----------------------------------------------------------------------------------------------------
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
