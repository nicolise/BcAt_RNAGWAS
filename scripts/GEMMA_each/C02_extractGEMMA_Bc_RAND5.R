#Nicole E Soltis
#07/02/18

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#Bc transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/04_GEMMAout/
#kmat script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A03_GEMMA_kmat.sh
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A04_runGEMMA_kmat_3genos.sh

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/")


#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")
mytime <- Sys.time()
  #each phenotype

for (k in 2:5){
  for (i in 1:9267){
    #actually: 1:9267 for Bc
    Sys.time()
    my_gemma <- read.table(paste("04_GEMMAout/col0_5Rand/rand",k,"/col0_MAF20NA10_RAND0",k,"_",i,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 4 seconds to read 1 phenotype
    #times 24000 = 27 hours
    #take top 10 SNP/phenotype
    #also save top 1 SNP/ phenotype
    my_gemma$pheno <- i
    my_gemma_top100 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top100 <- my_gemma_top100[1:100,]
    row.names(my_gemma_top100) <- c(1:100) + (100*(i-1))
    my_gemma_top10 <- my_gemma_top100[1:10,]
    my_gemma_top1 <- my_gemma_top10[1,]
    #and z scaling
    my_gemma.z <- my_gemma
    my_gemma.z$beta_z <- scale(my_gemma.z$beta, center = TRUE, scale = TRUE)
    my_gemma.z <- my_gemma.z[abs(my_gemma.z$beta_z) > 4,]
    mylgsnp <- my_gemma[abs(my_gemma$beta) > 0.5,]
    #this gives an error but it's fine
    try(ifelse( i == 1, write.table(my_gemma_top100, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top100SNPsample.txt"), sep = ",", col.names = TRUE), write.table(my_gemma_top100, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top100SNPsample.txt"), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top10, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top10SNPsample.txt"), sep = ",", col.names = TRUE), write.table(my_gemma_top10, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top10SNPsample.txt"), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top1, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top1SNPsample.txt"), sep = ",", col.names = TRUE), write.table(my_gemma_top1, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top1SNPsample.txt"), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma.z, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_topSNPsample_zscale.txt"), sep = ",", col.names = TRUE), write.table(my_gemma.z, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_topSNPsample_zscale.txt"), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(mylgsnp, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_top100_beta05SNP.txt"), sep = ",", col.names = TRUE), write.table(mylgsnp, paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",k,"_topSNP_beta05.txt"), sep = ",", col.names = FALSE, append = TRUE)))
    Sys.time()
  }
}
#}
mytime
Sys.time()