
#Nicole E Soltis 
#03/28/18
#add Phenos to FAM file for GEMMA
#----------------------------------------------------------------------
rm(list=ls())

#advice on permutations here:
#https://github.com/genetics-statistics/GEMMA/issues/93

#pipeline note:
#1. run C01_AtGWAS_PrepPhenos_lsmeans.R
#2. run C02_genostoPEDnMAP.R
#3. make sure there is a copy of plink executable in B05_GEMMA
#4. in command prompt: cd Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA/
#5. RUN ./plink --noweb --file C02_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE. NEXT STEP is customized by ogphenos/ permutation
#6. run this script (C03_prepPhenos.R)
#7. cd Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA/
#8. copy edited .fam, original .bim, .bed to C03_runGEMMA/
#9. copy bash script: cp scripts/GEMMA_lesions/norand_GEMMA_kmatrix.sh data/B05_GEMMA_les/
#10. cd to data/B05_GEMMA_les/
#11. calculate k-matrix with: bash norand_GEMMA_kmatrix.sh, mv files to D_03_kmat
#12. run GEMMA: bash norand_GEMMA_kmatrix_run.sh
#13. pheno order can be found in names(Phenos)

setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
Phenos <- read.csv("data/B05_GEMMA/C01_AtPrepFiles/lsmeans_zscale_allreads_At.csv")

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/B05_GEMMA/C02_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#col2 = V2 = Isolate
names(Phenos)[2] <- "V2"

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)


levels(Phenos_match$V2) <- c(levels(Phenos_match$V2), "MEAPGG")

#change Phenos_match$V2 "MEAP6G" to "MEAPGG"
Phenos_match[83, 1] <- "MEAPGG"

#remove 1.02.13, 94.4
##check these
Phenos_match <- Phenos_match[c(2:13,15:98),]

#double check now and save original order for FAM
myFAM_match$delete <- c(1:96)
setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#match and return to original order
myFAM_match2 <- merge(myFAM_match, Phenos_match, by="V2")
myFAM_match2 <- myFAM_match2[order(myFAM_match2$delete),]
#remove column "delete" and reorder column V1, V2
myFAM_match2 <- myFAM_match2[,c(2,1,3:6,8:ncol(myFAM_match2))]

#remove dummy phenotype column
myFAM_match2 <- myFAM_match2[,-c(6)]

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
write.csv(myFAM_match2, "B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_fam.csv")
write.table(myFAM_match2, "B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)
#if randomizing, copy _allphenos.fam to _rand.fam and also copy 01/.bed and 01/.bim to 02/_rand.
