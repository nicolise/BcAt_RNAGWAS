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