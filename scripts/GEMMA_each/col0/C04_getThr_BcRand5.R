#Nicole E Soltis
#11/05/18

#-------------------------------------------------------------------------------
#work on thresholding -- summarize across 5 permutations
#collect max value across 5 permuts for each SNP

#from GEMMA_lsm/D_GWAplots

rm(list=ls())
#start with a small file
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
#for now, don't care about annotating which phenotype was which (gene names for transcripts)
mydat_r1 <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_RAND1_top1SNPsample.txt")
mydat_r1 <- mydat_r1[,-c(2)]
mydat_r1$randrun <- 1
for (i in c(2:5)){
  mydat01 <- read.csv(paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=""))
  mydat01 <- mydat01[,-c(2)]
  mydat01$randrun <- i
  combdat <- mydat01[1,]
  for (j in c(1:9267)){
      if(mydat01[j,12] < mydat_r1[j,12]) {
        combdat[j,] <- mydat01[j,]
      } else { combdat[j,] <- mydat_r1[j,]}
  }
} 
  
  
#}


mydat <- mydat01
mydat$chr.snp <- paste(mydat$chr, mydat$ps, sep=".")

#save min p value per permutation
minrow <- mydat[mydat$p_score == min(mydat$p_score),]
#minrowall <- minrow
#do for all 5 permutations
minrowall <- rbind(minrowall, minrow)
minrowall$neg10p <- log10(minrowall$p_score)*-1
#highest = 7.53