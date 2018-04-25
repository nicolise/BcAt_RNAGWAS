
#Nicole E Soltis 
#03/28/18
#randomize phenotypes once each to try out permutations for GEMMA

#------------------------------------------------------------------------------
rm(list=ls())

#advice on permutations here:
#https://github.com/genetics-statistics/GEMMA/issues/93

#pipeline note:
#1. run D01_TABtoPEDnMAP.R
#2. copy plink executable to GEMMA_files
#3. in command prompt: cd to GEMMA_files
#4. RUN ./plink --noweb --file D_01_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE. NEXT STEP is customized by ogphenos/ permutation
#5. run this script (02_prepPhenos.R)
#6. cd to GEMMA_files
#7. copy edited _rand.fam, .bim, .bam to D_02_randGEMMA/
#8. calculate k-matrix with: bash fullrand_GEMMA_kmatrix.sh
#9. run GEMMA: bash fullrand_GEMMA_kmatrix_run.sh

#10: Pheno order is "

setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS")
Phenos <- read.csv("LSMeanCamLes4Map_FIN.csv")

setwd("~/Projects/BcSolGWAS/data")
setwd("~/Documents/GitRepos/BcSolGWAS/data/")
Phenos <- read.csv("GWAS_files/02_csvPrep/phenos/NewModel0711/BcSl_lsmeans_forbigRR.csv")
D.Phenos <- read.csv("GWAS_files/02_csvPrep/phenos/Domestication/BcSl_lsmeans_domest_forbigRR.csv")
#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("GEMMA_files/D_01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#first merge Phenos, D.Phenos
setdiff(Phenos$Igeno, D.Phenos$Igeno) #cool
myPhenos <- cbind(Phenos,D.Phenos)
myPhenos <- myPhenos[,-c(14)]
Phenos <- myPhenos
#col2 = V2 = Isolate
names(Phenos)[1] <- "V2"

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#add an empty variables 1.01.12 and 1.02.19 and Triple7 to Phenos
for (y in 2:length(Phenos_match)){
  Phenos_match[96,y] <- NA
  Phenos_match[97,y] <- NA
  Phenos_match[98,y] <- NA
}
levels(Phenos_match$V2) <- c(levels(Phenos_match$V2), "1.01.12", "1.02.19","Triple7", "MEAPGG")
Phenos_match[96,1] <- "1.01.12"
Phenos_match[97,1] <- "1.02.19"
Phenos_match[98,1] <- "Triple7"

#change Phenos_match$V2 "MEAP6G" to "MEAPGG"
Phenos_match[80, 1] <- "MEAPGG"

#remove 1.02.05, 94.4
Phenos_match <- Phenos_match[c(1:10,12:94,96:98),]

#reorder phenos
Phenos_match$V2
Phenos_match <- Phenos_match[c(1:5,94,6:15,95,16:44,93,45:91,96,92),]

#double check now and match orders
myFAM_match$delete <- c(1:96)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]
setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)
myFAM_match <- myFAM_match[1:96,]

#reorder phenotypes to make life way easier
names(Phenos_match)[13] <- "LA0480"
names(Phenos_match)[10] <- "LA0410"
Phenos_match <- Phenos_match[,c(1,10,13,2:9,11:12,14:16)]

myFAM_match2 <- cbind(myFAM_match, Phenos_match)
myFAM_match2 <- myFAM_match2[order(myFAM_match2$delete),]
myFAM_match2 <- myFAM_match2[,c(1:5,9:ncol(myFAM_match2))]

write.csv(myFAM_match2, "GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")
write.table(myFAM_match2, "GEMMA_files/D_02_randGEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)
#copy _allphenos.fam to _rand.fam and also copy 01/.bed and 01/.bim to 02/_rand.
