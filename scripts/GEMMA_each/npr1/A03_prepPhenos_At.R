
#Nicole E Soltis 
#03/28/18
#02_prepPhenos_At.R
#add Phenos to FAM file for GEMMA
#----------------------------------------------------------------------
rm(list=ls())

#pipeline note:
#1. copy files from GEMMA_eachAt_Bc/01_PLINK

setwd("/media/nesoltis/Soltis_npr1coi1_eQTL")
IsoNames <- read.csv("data/Wei_Bc/IsolateKey_Wei.csv")
MyReads <- read.csv("data/Wei_At/data_from_Wei_20160321/S2_ModelLsMeanAdjusted.csv")
MyReads <- MyReads[,-c(1)]

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/AtTranscript/01_PLINK/binMAF20NA10.fam")

#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.
## start here

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(23960,2,4:23959)]
names(MyReads)[1] <- "Isolate"

#first, split MyReads by plant accession. Then generate matching FAM files for each.
myread_coi1 <- MyReads[MyReads$HostGenotype=="coi.1",]
myread_col0 <- MyReads[MyReads$HostGenotype=="col.0",]
myread_npr1 <- MyReads[MyReads$HostGenotype=="npr.1",]

##do this for each of above
Phenos <- myread_npr1
#col2 = V2 = Isolate
Phenos <- Phenos[,c(1,3:length(Phenos))]
names(Phenos)[1] <- "V2"

Phenos_match <- Phenos[ order(Phenos$V2), ]
#remove non-genotyped 01.02.13 from Phenos
Phenos_match <- Phenos_match[!Phenos_match$V2 =="1.02.13",]
#and 94.4 and Control
Phenos_match <- Phenos_match[!Phenos_match$V2 =="94.4",]
Phenos_match <- Phenos_match[!Phenos_match$V2 =="Control",]
myFAM_match <- myFAM
myFAM_match$delete <- c(1:95)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]

## check that these are 0, 0, 95
setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#now add Phenos_match onto myFAM_match
myFAM_match2 <- merge(myFAM_match, Phenos_match, by="V2")
myFAM_match2 <- myFAM_match2[order(myFAM_match2$delete),]
#remove dummy phenotype (column 6)
#and reorder V1:V5
myFAM_match2 <- myFAM_match2[,c(2,1,3:5,8:length(myFAM_match2))]

##be sure to move new *.fam into the correct directory!
Sys.time()
write.table(myFAM_match2, "data/AtTranscript/npr1/02_GEMMA/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
Sys.time()