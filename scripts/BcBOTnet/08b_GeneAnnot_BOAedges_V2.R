#Nicole E Soltis
#07/19/18

#--------------------------------------------------------------------------
#this is for the BcBOTnet genes (Boa, Bot, Net5, could add ABA later)

rm(list=ls())
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
my.gtf <- my.gtf[,c(1:8,10,12,14)]
my.gtf.c1 <- my.gtf[my.gtf$V1=="Chromosome1",]

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv")
#-----------------------------------------------------------------------
#BOA gene annotation
my.boa <- mySNP_boa_named 
(my.boa[3,2]) #first SNP
(my.boa[3,length(my.boa)]) #last SNP
blah <- my.gtf.c1[my.gtf.c1$V4 < 5133,]
max(blah$V4) #none, so start at first SNP:
##so, beginning of Chr1 gene annotation is BOA cluster. Examine SNPs up to 2 kb before BOA SNP start = pos 5133 - 2e3
min(my.gtf.c1$V4)
blah <- my.gtf.c1[my.gtf.c1$V5 > 61757,]
View(blah)
##BcBoa17 goes up to 69449. 3 genes downstream is gene:Bcin01g00190. end position is 97414. 
#so, I'll select all SNPs in this region for excel plot to find deletion endpoints

#-------------------------------------------------------------------------------------
#select needed SNPs for phylogeny here
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
chr1.ped <- read.table("02b_Haploview/binMAF20NA10_chr1.ped")
row.names(chr1.ped) <- chr1.ped$V2
#index should link PED to SNP positions -- but each SNP twice in PED
chr1.info <- read.table("02b_Haploview/binMAF20NA10_chr1.info")

#approx region BoA cluster with surrounding edges: 3133 - 97414
mySNP_boa <- chr1.ped #lowest position is 4029- start here
myInfo_boa <- chr1.info

#duplicate SNPs in myInfo for matching with mySNP
myInfo_boa_db <- myInfo_boa[rep(row.names(myInfo_boa), 2), 1:2]
myInfo_boa_db <- myInfo_boa_db[order(myInfo_boa_db$V2),]
myInfo_boa_db$SNPcount <- c(1:nrow(myInfo_boa_db))
myInfo_boa_db$SNPnum <- paste("V",(myInfo_boa_db$SNPcount+6), sep="")
#index and PED matched
boa_names <- mySNP_boa[,1:6]
boa_SNPList <- myInfo_boa_db[,c(4,2)]
boa_SNPList <- as.data.frame(t(boa_SNPList))
for (i in 1:length(boa_SNPList)){
  levels(boa_SNPList[,i])<- c(levels(boa_SNPList[,i]),0,1,2)
}
colnames(boa_SNPList) <- as.character(unlist(boa_SNPList[1,]))
mySNP_boa_named <- mySNP_boa[,-c(1:6)]
mySNP_boa_named <- rbind(boa_SNPList, mySNP_boa_named)

#approx region BoA cluster with surrounding edges: 3133 - 97414
#lowest position is 4029- start here
boa_named_crop <- mySNP_boa_named[,1:500]
boa_named_crop2 <- as.data.frame(t(boa_named_crop))
boa_named_crop2$V2 <- as.integer(paste(boa_named_crop2$V2))
#next SNP after 97414 is 97730
boa_named_crop2 <- boa_named_crop2[boa_named_crop2$V2 <= 97730,]
boa_named_crop2 <- as.data.frame(t(boa_named_crop2))

boa_named_crop[,length(boa_named_crop)]
## this is a binary file from .ped. Includes duplicates of each SNP to fake diploidy. Includes SNP names (locations) from .info. Does not have genes mapped onto SNP locations.
write.csv(boa_named_crop, "02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa.csv")
#now add in gene locations!!

#---------------------------------------------------------------------------
my.gtf.boa <- my.gtf.c1[my.gtf.c1$V5 < 97730,]
##BcBoa17 goes up to 69449. 3 genes downstream is gene:Bcin01g00190. end position is 97414. 
my.genes.boa <- my.gtf.boa[,c("V4","V5","V12")]
#get gene start and stop for each
library("plyr")
#use this to hand-annotate the genes onto the excel plot
gene.ends.boa <- ddply(my.genes.boa, c("V12"), summarise,
                       geneMin = min(V4),
                       geneMax = max(V5))

#-----------------------------------------------------------------
#from excel plot:
#last base of deletion = 82614
#first base of deletion = first base of C1 = 4029
#now: pvclust with these new endpoints and vioplot with the new pvclust tree.
library(dplyr)

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_fulldel <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa.csv")
mySNP_boa_fulldel_gene <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa_wgenes.csv")
#remove SNPs outside of genes
#clustering includes only genic regions!
boaSNP_pvclust_geneonly <- mySNP_boa_fulldel_gene %>% 
  select_if(~ !any(is.na(.)))

##check which of these I am using for downstream analyses
#boaSNP_pvclust <- mySNP_boa_fulldel
boaSNP_pvclust <- boaSNP_pvclust_geneonly

boaSNP_pvclust <- as.data.frame(t(boaSNP_pvclust))
colnames(boaSNP_pvclust) = as.character(unlist(boaSNP_pvclust[1,]))
boaSNP_pvclust <- boaSNP_pvclust[-c(1),]
boaSNP_pvclust$V2 <- as.integer(paste(boaSNP_pvclust$V2))
# #This step does not make sense- should cluster including full genic region!!
# #crop to only include deletion
# boaSNP_pvclust <- boaSNP_pvclust[boaSNP_pvclust$V2<=82614,]
#get rid of SNPnum and V2 for pvclust
boaSNP_pvclust <- boaSNP_pvclust[,-c(1:2)]

library(pvclust)
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

#now run pvclust
mytime <- Sys.time()
pvclust_boa_out <- pvclust(boaSNP_myvar, method.dist = "cor", method.hclust= "average", nboot=1000)
mytime
Sys.time()

#with PV method 
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/PVclusters/BotcynicAcid_IsolatesPVClust_fulldel_onlygenes.jpg", sep=""), width=30, height=10, units='in', res=600)
plot(pvclust_boa_out)
pvrect(pvclust_boa_out, alpha=0.95)
dev.off()

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
BoaClusts <- myclusters[,c("Isolate","boa_pv_fulldel")] #can also do with boa_pv_fulldel_onlygene
BoaPhenosPV <- merge(BoaPhenos,BoaClusts, by = "Isolate")

#plot mean phenotype by cluster
library(ggplot2)
p <- ggplot(BoaPhenosPV, aes(factor(boa_pv_fulldel), mean.Pheno))
setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/Cluster_Means/BotcynicAcid_PVclust_fulldel.jpg", sep=""), width=8, height=6, units='in', res=600)
p + geom_violin(trim=FALSE, draw_quantiles = c(0.25, 0.5, 0.75), fill="slateblue1", alpha=0.4) + geom_jitter(height = 0, width = 0.1) + 
  geom_text(aes(label=Isolate), hjust=-1.5, vjust=0, size=2)+ 
  labs(x="Cluster Membership", y = "Mean Expression Across Network")+ 
  theme_bw()
dev.off()
