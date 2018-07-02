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

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA")
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

#bot limits C12 2,217,306	2,245,431
#min
myInfoc12.b <- myInfoc12[myInfoc12$V4 < 2245431,]
mymax <- max(myInfoc12.b$V4)
myminrow <- which(grepl(mymax, myInfoc12$V4))
#calculate which column for each in ped
mymincol <- 2*myminrow + 6
#max
myInfoc12.b <- myInfoc12[myInfoc12$V4 > 2217306,]
mymin <- min(myInfoc12.b$V4)
mymaxrow <- which(grepl(mymin, myInfoc12$V4))
#calculate which column for each in ped
mymaxcol <- 2*mymaxrow + 6
#keep column: 1:6, 28702:28228
binPEDc12_bot <- binPEDc12[,c(1:6, 28228:28702)]
myInfoc12_bot <- myInfoc12[c(14111:14348),]
#remove last 2 snps from each for haploview - appear to be outside high LD region
#rm snp 2243688
binPEDc12_bot <- binPEDc12[,c(1:6, 28228:28699)]
myInfoc12_bot <- myInfoc12[c(14111:14345),]

write.table(binPEDc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_crop.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_crop.info", col.names=FALSE, row.names=FALSE)


binPEDc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_crop.ped")
myInfoc12_bot <- read.table("02b_Haploview/binMAF20NA10_chr12_bot_crop.info")
row.names(binPEDc12_bot) <- binPEDc12_bot$V2
clusters_bot <- hclust(dist(binPEDc12_bot[, 7:478]))
plot(clusters_bot)
clusterCut_bot <- cutree(clusters_bot, 3)
#cluster membership of each isolate
table(clusterCut_bot, binPEDc12_bot$V2)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/HClust/Botrydial_IsolatesCluster.jpg", sep=""), width=8, height=6, units='in', res=600)
plot(clusters_bot, cex=0.5)
dev.off()

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

