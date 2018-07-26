#Nicole E Soltis
#02/21/18
#02_FAMaddPhenos

#----------------------------------------------------------------------------
rm(list=ls())

#pipeline note:
#1. run A01_TABtoPEDnMAP_rmIso.R
#2. make sure there is a copy of plink executable in data/GEMMA_eachAt_Bc
#3. in command prompt: cd Documents/GitRepos/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc
#4. RUN ./plink --noweb --file 01_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } do this ONCE. NEXT STEP is customized by ogphenos/ permutation
#5. copy these files to GEMMA_eachBc_At/01_PLINK
#6. run this script (A02_prepPhenos_Bc.R)

##start here
#7. cd Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA/
#8. copy edited .fam, original .bim, .bed to C03_runGEMMA/
#9. copy bash script: cp scripts/GEMMA_lesions/norand_GEMMA_kmatrix.sh data/B05_GEMMA_les/
#10. cd to data/B05_GEMMA_les/
#11. calculate k-matrix with: bash C04_runGEMMA_allAt_kmat.sh, mv files to C04_kmat
#12. move the whole thing to Data/ drive
#13. on Data/ run GEMMA: bash C05_runGEMMA_allAt_kmat_run.sh
#14. pheno order can be found in names(Phenos)

setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
IsoNames <- read.csv("data/Vivian_Bc/IsolateKey_Vivian.csv")
MyReads <- read.csv("data/Vivian_Bc/result.lsm.csv") #9270 phenotypes
MyReads <- MyReads[,-c(1)]
#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(9270,2:9269)]
names(MyReads)[1] <- "Isolate"

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/GEMMA_eachAt_Bc/01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

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
write.table(myFAM_match2, "data/GEMMA_eachAt_Bc/02_GEMMA/binMAF20NA10.fam", row.names=FALSE, col.names=TRUE)
Sys.time()

myFAM_check <- read.table("data/GEMMA_eachAt_Bc/col0/02_GEMMA/binMAF20NA10.fam")
myFAM_check2 <- read.table("data/GEMMA_eachBc_At/col0/02_GEMMA/binMAF20NA10.fam")
