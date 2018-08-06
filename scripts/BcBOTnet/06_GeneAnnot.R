#Nicole E Soltis
#07/19/18

#--------------------------------------------------------------------------
#this is for the BcBOTnet genes (Boa, Bot, Net5, could add ABA later)

rm(list=ls())
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)

#snp location:
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
#bot limits C12 2,217,306	2,245,431
## where does this come from?
#index should link PED to SNP positions -- but each SNP twice in PED
binPEDc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_recrop.ped")
row.names(binPEDc12_bot) <- binPEDc12_bot$V2
binPEDc12_info <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_recrop.info")
#approx region Bot deletion : C12 SNP 240 - 440
mySNP_bot <- binPEDc12_bot
head(binPEDc12_info)
myInfo_bot <- binPEDc12_info
myInfo_db <- myInfo_bot[rep(row.names(myInfo_bot), 2), 1:2]
myInfo_db <- myInfo_db[order(myInfo_db$V2),]
myInfo_db$SNPcount <- c(1:nrow(myInfo_db))
myInfo_db$SNPnum <- paste("V",(myInfo_db$SNPcount+6), sep="")
#index and PED matched
bot_names <- mySNP_bot[,1:6]
bot_SNPList <- myInfo_db[,c(4,2)]
bot_SNPList <- as.data.frame(t(bot_SNPList))
for (i in 1:length(bot_SNPList)){
  levels(bot_SNPList[,i])<- c(levels(bot_SNPList[,i]),0,1,2)
}
colnames(bot_SNPList) <- as.character(unlist(bot_SNPList[1,]))
mySNP_bot_named <- mySNP_bot[,-c(1:6)]
mySNP_bot_named <- rbind(bot_SNPList, mySNP_bot_named)

#this now has an added row of gene id- copied by hand from Excel .xls version
#write.csv(mySNP_bot_named, "02b_Haploview/binMAF20NA10_chr12_bot_recrop_named.csv")
#now add in gene locations!!
mySNP_bot_named <- read.csv("02b_Haploview/binMAF20NA10_chr12_bot_recrop_named.csv")
#-----------------------------------------------------------------------
#BOT gene annotation
num.genes <- my.gtf[unique(my.gtf$V12),]

my.gtf <- my.gtf[,c(1:8,10,12,14)]

my.gtf.c1 <- my.gtf[my.gtf$V1=="Chromosome1",]
my.gtf.c12 <- my.gtf[my.gtf$V1=="Chromosome12",]

#bot gene matching
my.bot <- mySNP_bot_named 
(my.bot[2,1]) #first SNP
(my.bot[2,length(my.bot)]) #last SNP
blah <- my.gtf.c12[my.gtf.c12$V4 < 2217400,]
max(blah$V4)
my.gtf.bot <- my.gtf.c12[my.gtf.c12$V4 > 2217305,]
blah <- my.gtf.c12[my.gtf.c12$V5 > 2243688,]
min(blah$V5)
my.gtf.bot <- my.gtf.bot[my.gtf.bot$V5 < 2244208,]

my.genes.bot <- my.gtf.bot[,c("V4","V5","V12")]
#get gene start and stop for each
library("plyr")
gene.ends.bot <- ddply(my.genes.bot, c("V12"), summarise,
                 geneMin = min(V4),
                geneMax = max(V5))
#-----------------------------------------------------------------------
#boa
binPEDc1_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa.ped")
row.names(binPEDc1_boa) <- binPEDc1_boa$V2

#snp location:
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
#boa limits
## where does this come from?
#index should link PED to SNP positions -- but each SNP twice in PED
binPEDc1_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa_fix.ped")
row.names(binPEDc1_boa) <- binPEDc1_boa$V2
info_boa <- read.table("02b_Haploview/binMAF20NA10_chr1_boa_fix.info")
#approx region Bot deletion : C12 SNP 240 - 440
mySNP_boa <- binPEDc1_boa
head(info_boa)
myInfo_boa <- info_boa
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

#write.csv(mySNP_boa_named, "02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv")
#now add in gene locations!!
mySNP_boa_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv")
#-----------------------------------------------------------------------
#BOA gene annotation
my.boa <- mySNP_boa_named 
(my.boa[2,1]) #first SNP
(my.boa[2,length(my.boa)]) #last SNP
blah <- my.gtf.c1[my.gtf.c1$V4 < 5133,]
max(blah$V4) #none, so start at first SNP:
min(my.gtf.c1$V4)
my.gtf.boa <- my.gtf.c1[my.gtf.c1$V4 > 5428,]
blah <- my.gtf.c1[my.gtf.c1$V5 > 61757,]
min(blah$V5) 
my.gtf.boa <- my.gtf.boa[my.gtf.boa$V5 < 63780,]

my.genes.boa <- my.gtf.boa[,c("V4","V5","V12")]
#get gene start and stop for each
library("plyr")
gene.ends.boa <- ddply(my.genes.boa, c("V12"), summarise,
                       geneMin = min(V4),
                       geneMax = max(V5))
#---------------------------------------------------------------------------
#net5
binPEDc1_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5.ped")
row.names(binPEDc1_net5) <- binPEDc1_net5$V2

#snp location:
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
#net5 limits
#index should link PED to SNP positions -- but each SNP twice in PED
binPEDc1_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5_fix.ped")
row.names(binPEDc1_net5) <- binPEDc1_net5$V2
info_net5 <- read.table("02b_Haploview/binMAF20NA10_chr1_net5_fix.info")
#approx region Bot deletion : C12 SNP 240 - 440
mySNP_net5 <- binPEDc1_net5
head(info_net5)
myInfo_net5 <- info_net5
#duplicate SNPs in myInfo for matching with mySNP
myInfo_net5_db <- myInfo_net5[rep(row.names(myInfo_net5), 2), 1:2]
myInfo_net5_db <- myInfo_net5_db[order(myInfo_net5_db$V2),]
#name SNPs with V7...
myInfo_net5_db$SNPcount <- c(1:nrow(myInfo_net5_db))
myInfo_net5_db$SNPnum <- paste("V",(myInfo_net5_db$SNPcount+6), sep="")
#index and PED matched
net5_names <- mySNP_net5[,1:6]
net5_SNPList <- myInfo_net5_db[,c(4,2)]
net5_SNPList <- as.data.frame(t(net5_SNPList))
for (i in 1:length(net5_SNPList)){
  levels(net5_SNPList[,i])<- c(levels(net5_SNPList[,i]),0,1,2)
}
colnames(net5_SNPList) <- as.character(unlist(net5_SNPList[1,]))
mySNP_net5_named <- mySNP_net5[,-c(1:6)]
mySNP_net5_named <- rbind(net5_SNPList, mySNP_net5_named)

#write.csv(mySNP_net5_named, "02b_Haploview/binMAF20NA10_chr1_net5_fix_named.csv")
#now add in gene locations!!
mySNP_net5_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_net5_fix_named.csv")
#-----------------------------------------------------------------------
#NET5 gene annotation
my.net5 <- mySNP_net5_named 
(my.net5[2,1]) #first SNP
(my.net5[2,length(my.net5)]) #last SNP
blah <- my.gtf.c1[my.gtf.c1$V4 < 4026366,]
max(blah$V4) #4023303
my.gtf.net5 <- my.gtf.c1[my.gtf.c1$V4 > 4023303,]
blah <- my.gtf.c1[my.gtf.c1$V5 > 4073877,]
min(blah$V5) 
my.gtf.net5 <- my.gtf.net5[my.gtf.net5$V5 < 4107137,]

my.genes.net5 <- my.gtf.net5[,c("V4","V5","V12")]
#get gene start and stop for each
library("plyr")
gene.ends.net5 <- ddply(my.genes.net5, c("V12"), summarise,
                       geneMin = min(V4),
                       geneMax = max(V5))