
#Nicole E Soltis 
#03/28/18
#randomize phenotypes once each to try out permutations for GEMMA

#----------------------------------------------------------------------
rm(list=ls())

#advice on permutations here:
#https://github.com/genetics-statistics/GEMMA/issues/93

#pipeline note:
#1. run D01_TABtoPEDnMAP.R
#2. copy plink executable to B05_GEMMA_les
#3. in command prompt: cd to B05_GEMMA_les
#4. RUN ./plink --noweb --file D_01_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE. NEXT STEP is customized by ogphenos/ permutation
#5. run this script (02_prepPhenos.R)
#6. cd to GEMMA_files
#7. copy edited .fam, original .bim, .bed to D_02_randGEMMA/
#8. copy bash script: cp scripts/GEMMA_lesions/norand_GEMMA_kmatrix.sh data/B05_GEMMA_les/
#9. cd to data/B05_GEMMA_les/
#9. calculate k-matrix with: bash norand_GEMMA_kmatrix.sh, mv files to D_03_kmat
#10. run GEMMA: bash norand_GEMMA_kmatrix_run.sh
#11: Pheno order to run is: "Col0.Les","coi1.Les","npr1.Les","pad3.Les","tga3.Les","anac055.Les"

setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
setwd("~/Documents/GitRepos/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
Phenos <- read.csv("LSMeanCamLes4Map_FIN.csv")

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/")
#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("B05_GEMMA_les/D_01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#col2 = V2 = Isolate
names(Phenos)[2] <- "V2"
#keep only lesion phenos since that's what we want to run
Phenos <- Phenos[,c("V2","Col0.Les","coi1.Les","npr1.Les","pad3.Les","tga3.Les","anac055.Les")]

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#add an empty variable 1.01.12 to Phenos
for (y in 2:length(Phenos_match)){
  Phenos_match[98,y] <- NA
}

levels(Phenos_match$V2) <- c(levels(Phenos_match$V2), "1.01.12", "MEAPGG")
Phenos_match[98,1] <- "1.01.12"

#change Phenos_match$V2 "MEAP6G" to "MEAPGG"
Phenos_match[83, 1] <- "MEAPGG"

#remove 1.02.13, 94.4
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

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
write.csv(myFAM_match2, "B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_fam.csv")
write.table(myFAM_match2, "B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_allphenos.fam", row.names=FALSE, col.names=TRUE)
#if randomizing, copy _allphenos.fam to _rand.fam and also copy 01/.bed and 01/.bim to 02/_rand.
