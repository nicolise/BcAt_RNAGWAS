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
#binMAPc1 is 23350 SNPs which should be the first (23350 * 2) + 6 = 46706 columns of binPED (binPED V6 is NULL PHENOTYPE)
binPEDc1 <- binPED[,c(1:46706)]
write.table(binPEDc1, "02b_Haploview/binMAF20NA10_chr1.ped", col.names=FALSE, row.names=FALSE)
myInfoc1 <- binMAPc1[,c(2,4)]
write.table(myInfoc1, "02b_Haploview/binMAF20NA10_chr1.info", col.names=FALSE, row.names=FALSE)

binPEDc1 <- read.table("02b_Haploview/binMAF20NA10_chr1.ped")
myInfoc1 <- read.table("02b_Haploview/binMAF20NA10_chr1.info")

#boa limits C1 5,429	61,522
##where does this come from?
#find which row for each limit in info
#min, START of boa cluster
myInfoc1.b <- myInfoc1[myInfoc1$V2 < 5429,]
#closest SNP outside region is 5133
#this only works if exact match
max(myInfoc1.b$V2) #nearest SNP BELOW 5429
which(grepl(5133, myInfoc1$V2)) #row 11 is START for boa cluster
myInfoc1.b[myInfoc1.b$V2==5133,] #double check, row 11
#calculate which column for each in ped
2*11 + 6
#keep PED columns: 1:6, 28:x
#max, STOP of boa cluster
myInfoc1.b <- myInfoc1[myInfoc1$V2 > 61522,]
#closest SNP outside region is 61757
#this only works if exact match
min(myInfoc1.b$V2) #nearest SNP ABOVE 61522
which(grepl(61757, myInfoc1$V2)) #row 133 is STOP for boa cluster
myInfoc1.b[myInfoc1.b$V2==61757,] #double check, row 133
#calculate which column for each in ped
2*133 + 6
#keep PED columns: 1:6, 28:272
#28 is SECOND column of SNP 11-- need from 27
#in fix file, starts 27. in old file, starts 28. 
binPEDc1_boa <- binPEDc1[,c(1:6, 27:272)] #select PED columns
myInfoc1_boa <- myInfoc1[11:133,] #select Info rows
#check myInfo vs. binPED
123*2+6 
write.table(binPEDc1_boa, "02b_Haploview/binMAF20NA10_chr1_boa_fix.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_boa, "02b_Haploview/binMAF20NA10_chr1_boa_fix.info", col.names=FALSE, row.names=FALSE)

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
#start of first SNP: for net5 is 46381! old file goes from 46382
binPEDc1_net5 <- binPEDc1[,c(1:6, 46381:46604)]
myInfoc1_net5 <- myInfoc1[23188:23299,]
#check length binPED
112*2+6
write.table(binPEDc1_net5, "02b_Haploview/binMAF20NA10_chr1_net5_fix.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc1_net5, "02b_Haploview/binMAF20NA10_chr1_net5_fix.info", col.names=FALSE, row.names=FALSE)

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

binPEDc12 <- read.table("02b_Haploview/binMAF20NA10_chr12.ped") #30083
myInfoc12 <- read.table("02b_Haploview/binMAF20NA10_chr12.info") #15039 * 2 = 30078

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
binPEDc12_bot <- binPEDc12[,c(1:6, 28228:28699)] #first and last SNPs have pair- good!
myInfoc12_bot <- myInfoc12[c(14111:14346),]
#check size binPED
236*2+6
##size is correct for recrop, PED slightly off for crop.ped
write.table(binPEDc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_recrop.ped", col.names=FALSE, row.names=FALSE)
write.table(myInfoc12_bot, "02b_Haploview/binMAF20NA10_chr12_bot_recrop.info", col.names=FALSE, row.names=FALSE)
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
