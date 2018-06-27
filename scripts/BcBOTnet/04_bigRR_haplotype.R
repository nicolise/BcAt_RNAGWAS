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
#try just C1 and C12 for BotBoaNet5
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA")
myPED <- read.table("01_PLINK/dpcharMAF20NA10.ped")
myMAP <- read.table("01_PLINK/dpcharMAF20NA10.map")

#try haploview format:
#*.info is 2 tab-separated columns with no headers. C1 is SNP names, C2 is SNP positions
#*.ped is a normal .ped file with no headers. 
#6 columns of labels, then paired SNP states

myInfo <- myMAP[,c(2,4)]
write.table(myInfo, "02b_Haploview/charMAF20NA10.info", col.names=FALSE, row.names=FALSE)

write.table(myPED, "02b_Haploview/charMAF20NA10.ped", col.names=FALSE, row.names=FALSE)

#haploview error: >2 alleles at marker 1987 -- use binary coded data!
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

